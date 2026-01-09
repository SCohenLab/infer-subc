from typing import List, Union
from pathlib import Path
import warnings
import time

import numpy as np
import pandas as pd

from infer_subc.quantification.morphology import get_morphology_metrics
from infer_subc.quantification.batch import load_existing_keys_csv, append_atomic_csv
from infer_subc.utils.batch import list_image_files, find_segmentation_tiff_files
from infer_subc.core.file_io import read_czi_image, read_tiff_image


def get_regions_morphology(source_file_path: str,
                           list_region_names: Union[List[str], None]=None,
                           list_region_segs: Union[List[np.ndarray], None]=None,
                           list_intensity_img: Union[List[np.ndarray], None]=None,
                           list_channel_names: Union[List[str], None]=None,
                           mask_name: Union[str, None]=None,
                           scale: Union[tuple, None]=None) -> pd.DataFrame:
    """
    Measure morphology metrics of masks/regions included in the large infer-subc pipeline (e.g. cell, nucleus, etc.).

    Parameters
    ------------
    source_file: str
        Path to the source image file. This will be used as part of the metadata information in the output table. 
        The input images are not derived from this path, but rather are provided directly as arrays in the list_obj_segs and 
        list_intensity_img variables below.
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
    list_intensity_img: Union[List[np.ndarray], None]
        List of 3D intensity channels from the raw image. Any number of channels can be included.
        These names will be used to rename the intensity measurement columns in the output table.
        If no intensity analysis is to be included, specify None here.
    list_channel_names: Union[List[str], None]
        List of names for each intensity channel provided in list_intensity_img. The order should match the order of the channels in list_intensity_img.
    mask_name: Union[str, None]
        Name of the region to use as the mask for analysis; if not specified, the entire image will be quantified.
        The mask_name should match one of the names provided in list_region_names. This object will be used to mask
        all other objects before quantitative analysis is performed. It will also be included as one of the analyzed objects.
    scale: Union[tuple,None] = None
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)
            
    Returns
    -------------
    pandas dataframe of containing regionprops measurements (columns) for each object in the segmentation image (rows) and the regionprops object

    """
    # Validate inputs
    if list_region_names is None or list_region_segs is None:
        raise ValueError("You must provide both list_region_names and list_region_segs arguments.")
    if len(list_region_names) != len(list_region_segs):
        raise ValueError("The length of list_region_names must match the length of list_region_segs.")
    if list_intensity_img is None or list_channel_names is None:
        raise ValueError("You must provide both list_intensity_img and list_channel_names arguments.")
    if len(list_intensity_img) != len(list_channel_names):
        raise ValueError("The length of list_intensity_img must match the length of list_channel_names.")

    if isinstance(source_file_path, str): source_file_path = Path(source_file_path)
    print(f"Quantifying region morphology from {source_file_path.name}")


    # specify the mask image to use during quantification
    if mask_name is None:
        mask = None
        print("No mask name provided. No mask will be applied before analysis.")
    elif mask_name not in list_region_names:
        print(f"Mask '{mask_name}' not found in `list_region_names`:{list_region_names}. No mask will be applied before analysis.")
        mask = None
        mask_name = None
    else:
        mask = list_region_segs[list_region_names.index(mask_name)]
        print(f"Using '{mask_name}' as the mask for analysis.")

    # merge intensity images to create a single np.ndarray
    if list_intensity_img is None:
        intensity_img = None
        print("No intensity images provided. Morphology metrics that require intensity images will not be calculated.")
    else:
        intensity_img = np.stack(list_intensity_img, axis=0)

    # empty list to collect a morphology data for each organelle
    regions_tab = []

    # loop through the list of organelles and run the get_morphology_metrics function
    for j, target in enumerate(list_region_names):
        region_seg = list_region_segs[j]

        # verify only one object per mask image; if more than one, combine them into a single object
        ## TODO: update to multi-object analysis later
        unique_objs = np.unique(region_seg)
        unique_objs = unique_objs[unique_objs != 0]  # exclude background
        if len(unique_objs) > 1:
            warnings.warn(f"More than one object found in region segmentation '{list_region_names[j]}'. Combining all objects into a single object for analysis.")
            region_seg = (region_seg > 0).astype(int)
        else:
            region_seg = region_seg

        # run get_morphology_metrics function to output a table of measurements
        region_metrics = get_morphology_metrics(segmentation_img=region_seg, 
                                            seg_name=target,
                                            intensity_img=intensity_img, 
                                            intensity_ch_names=list_channel_names,
                                            channel_axis=0, # default to 0 because intensities are merged from list on axis 0
                                            mask=mask,
                                            mask_name=mask_name,
                                            scale=scale)
        
        # add table to list above
        regions_tab.append(region_metrics)

    # combine the lists for each organelle into one table
    final_region_tab = pd.concat(regions_tab, ignore_index=True)

    # add a new column to list the name of the image these data are derived from 
    final_region_tab.insert(loc=0,column='image_name',value=source_file_path.stem)

    return final_region_tab



def batch_process_regions_morph(dataset_name: str,
                                 raw_path: Union[Path,str], 
                                 seg_path: Union[Path,str],
                                 quant_path: Union[Path, str], 
                                 raw_file_type: str,
                                 region_names: Union[List[str], None],
                                 channel_axis: Union[int, None]=None,
                                 channel_names: Union[List[int], None]=None,
                                 mask_name: Union[str, None]=None,
                                 use_scale: bool=True,
                                 seg_suffix: Union[str, None]=None):
    """  
    batch process quantification of the regions morphology for a single dataset. This function is currently 
    optimized to process images from one file folder per image type (e.g., raw, segmentation) the output csv 
    files are saved to the indicated quant_path folder.

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
    region_names: Union[List[str], None]=None
        List of region names to analyze. Usually ['cell', 'nuc'] for cell mask and nucleus.
        If no regions are to be included, specify None here.
    channel_axis : int
        Axis corresponding to the channels in the image data
    channel_names: List[str]
        List of channel names associated to each channel in the raw image data; if you wish to exclude a particular channel
        from the intensity analysis, write None instead of the channel name.
    mask_name: Union[str, None]=None
        Name of the region to use for segmentation (if any). This name should be included in the regions_name variable.
        If None, the entire image will be quantified.
    use_scale: bool=True
        Whether to apply scaling to the quantitative data; scaled data will be in real world units (e.g., microns) rather than pixels/voxels
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

    if region_names is None:
        raise ValueError("No region names provided. Please provide at least one region name to analyze.")

    # check if any existing data is present in outfiles to skip already processed images
    unique_keys = ['dataset', 'image_name']

    regions_path = quant_path / f"{dataset_name}_regions_morphology_metrics.csv"
    existing_morpho_keys = load_existing_keys_csv(regions_path, unique_keys)

    # reading list of files from the raw path
    img_file_list = list_image_files(raw_path, raw_file_type)
    len_file_list = len(img_file_list)

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
            filez = find_segmentation_tiff_files(img_f, region_names, seg_path, seg_suffix)

            # read in raw file and metadata
            img_data, meta_dict = read_czi_image(filez["raw"])

            # create intensities from raw file as list baseed on channel_name list
            if channel_names is None:
                intensities = None
                print("No intensity channel information provided.")
            else:
                if channel_axis != 0:
                    img_data = np.moveaxis(img_data, channel_axis, 0)
                intensities = [img_data[i] for i, ch in enumerate(channel_names) if ch is not None]
                channel_names = [ch for ch in channel_names if ch is not None]

            # store region images as list
            regions = [read_tiff_image(filez[org]) for org in region_names]

            # define the scale
            if use_scale is True:
                scale = meta_dict['scale']
            else:
                scale = None

            regions_metrics = get_regions_morphology(source_file_path=img_f,
                                                  list_region_names=region_names,
                                                  list_region_segs=regions,
                                                  list_intensity_img=intensities,
                                                  list_channel_names=channel_names,
                                                  mask_name=mask_name,
                                                  scale=scale)
            
            # save the morphology table data per image directly to csv
            regions_metrics.insert(loc=0,column='dataset',value=dataset_name)
            append_atomic_csv(regions_path, regions_metrics)
            del regions_metrics  # free up memory

            end2 = time.time()
            print(f"Completed quantification of {meta_dict['file_name']} in {(end2-img_start)/60} mins.")
            print(f"{count}/{len_file_list} images have been processed.")
            print(f"Time elapsed: {(end2-img_start)/60} mins")

    end = time.time()
    print(f"Quantification for {count} files is COMPLETE! Files saved to '{quant_path}'.")
    print(f"It took {(end - start)/60} minutes to quantify these files.")



def batch_regions_morph_summary_stats(csv_path_list: List[str],
                                        out_path: str,
                                        out_prefix: str,
                                        region_names: List[str]):
    """" 
    csv_path_list: List[str],
        A list of path strings where .csv files to analyze are located.
    out_path: str,
        A path string where the summary data file will be output to
    out_prefix: str
        The prefix used to name the output file. An "_" will be included between this prefix and the file suffix.
    region_names: List[str],
        A list of region names used in the region morphology quantification (batch_process_region_morph function)
    """
    # for keeping track of dataset and file numbers
    ds_count = 0
    fl_count = 0

    ###################
    # Read in the csv files and combine them into one of each type
    ###################
    # create empty list to hold the regions tables from different experiments
    regions_tab = []

    # loop through all of the locations listed above and find the _regions_morph files; append them to the list above
    for loc in csv_path_list:
        # list all csv files in the location
        files_store = sorted(loc.glob("*.csv"))

        # find the unique datasets in this location based on the prefixes before "_regions_morphology_metrics"
        prefixes = set(f.name.split("_regions_morphology_metrics")[0] for f in files_store if "_regions_morphology_metrics" in f.name)
        for prefix in prefixes:
            ds_count += 1
            # select only the files from this dataset
            files_subset = [f for f in files_store if f.name.startswith(prefix +"_regions_morphology_metrics")]
            for file in files_subset:
                fl_count += 1
                stem = file.stem
                if "_regions_morph" in stem:
                    test_regions = pd.read_csv(file, index_col=0)
                    regions_tab.append(test_regions)

    # combine the regions_morph lists found above into one table
    regions_df = pd.concat(regions_tab,axis=0, join='outer').reset_index()
    print(f"Found {fl_count} files from {ds_count} dataset(s) across {len(csv_path_list)} location(s).")


    ###################
    # summary stat group
    ###################
    group_by = ['dataset', 'image_name', 'mask_name', 'scale', 'object']
    sharedcolumns = ["SA_to_volume_ratio", "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length"]  + list(regions_df.filter(regex=".*intensity.*").columns)
    ag_func_standard = ['mean', 'median', 'std']

    ###################
    # summarize morphology metrics
    ###################
    tab1 = regions_df[group_by + ['label']].groupby(group_by).agg(['count'])
    tab1.rename(columns={'label': 'region'}, inplace=True)
    tab2 = regions_df[group_by + ['volume', 'surface_area']].groupby(group_by).agg(['sum'] + ag_func_standard)
    tab3 = regions_df[group_by + sharedcolumns].groupby(group_by).agg(ag_func_standard)
    regions_summary = pd.merge(tab1, tab2, 'outer', on=group_by)
    regions_summary = pd.merge(regions_summary, tab3, 'outer', on=group_by)

    # Get mask_name and corresponding volume column per group & calculate volume fraction
    mask_names = regions_df.groupby(group_by)['mask_name'].first()
    mask_volume_data = regions_df.groupby(group_by).first().apply(lambda row: row[f"{mask_names.loc[row.name]}_volume"], axis=1)
    regions_summary.insert(regions_summary.columns.get_loc(('volume', 'sum')) + 1, ('volume', 'fraction'), regions_summary[('volume', 'sum')]/mask_volume_data)

    #######################
    # fill gaps & NA values
    #######################
    # Ensure all possible interactions are represented (if missing fill with NaN)
    for ind in regions_summary.index.droplevel(4).unique().to_list():
        for row in region_names:
            if ind+(row,) not in regions_summary.index:
                regions_summary.loc[ind+(row,)] = np.nan

    # fill NA with 0 for specific columns
    fill_dict = {('region', 'count'): 0, 
                ('volume', 'sum'): 0,
                ('surface_area', 'sum'): 0,
                ('volume', 'fraction'): 0}
    regions_summary = regions_summary.fillna(value=fill_dict)

    # if (region, count) is 1, set mean, median, and std to NaN
    single_site_mask = regions_summary[('region', 'count')] == 1
    for col in sharedcolumns+['volume', 'surface_area']:
        regions_summary.loc[single_site_mask, (col, 'std')] = np.nan

    regions_summary.sort_index(inplace=True)

    ###################
    # flatten datasheet and export
    ###################
    # export before unstacking
    if (Path(out_path) / f"{out_prefix}_per_region_morphology_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}_per_region_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
    else:
        regions_summary.to_csv(str(out_path) + f"/{out_prefix}_per_region_morphology_summarystats.csv", mode='x')
        print(f"Exported per-region morphology summary statistics (before unstacking) to {out_path}/{out_prefix}_per_region_morphology_summarystats.csv")
    regions_morph_final = regions_summary.unstack(-1)
    regions_morph_final.columns = ["_".join((col_name[1], col_name[-1], col_name[0])) for col_name in regions_morph_final.columns.to_flat_index()]
    regions_morph_final.columns = [col.replace('sum', 'total') for col in regions_morph_final.columns]
    regions_morph_final.columns = [col.replace(col, 'mask_volume') if 'mask' in col else col for col in regions_morph_final.columns]
    regions_morph_final = regions_morph_final.loc[:, ~regions_morph_final.columns.duplicated()]
    regions_morph_final.reset_index(inplace=True)

    ###################
    # export summary sheets
    ###################
    if (Path(out_path) / f"{out_prefix}_regions_morphology_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}_regions_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
    else:
        regions_morph_final.to_csv(str(out_path) + f"/{out_prefix}_regions_morphology_summarystats.csv", mode='x')
        print(f"Exported regions morphology summary statistics (after unstacking) to {out_path}/{out_prefix}_regions_morphology_summarystats.csv")
    print(f"Regions morphology summary is complete.")
    return regions_summary