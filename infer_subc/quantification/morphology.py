from typing import List, Union
from pathlib import Path
import warnings
import time

import numpy as np
import pandas as pd
from skimage.measure import regionprops_table

from infer_subc.core.img import *
from infer_subc.organelles import * 

from infer_subc.utils.batch import list_image_files, find_segmentation_tiff_files
from infer_subc.core.file_io import read_czi_image, read_tiff_image


# Apply skimage regionprops function to quantify the morphology of a single object type from a single cell
def get_morphology_metrics(segmentation_img: np.ndarray, 
                           seg_name: str, 
                           intensity_img, 
                           mask: np.ndarray, 
                           mask_name: str,
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
    mask:
        a binary np.ndarray mask of the area to measure from; this image should be the same shape as the segmentation file
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
    input_labels = _assert_uint16_labels(segmentation_img)
    input_labels = apply_mask(input_labels, mask)

    ##########################################
    ## CREATE LIST OF REGIONPROPS MEASUREMENTS
    ##########################################
    # start with LABEL
    properties = ["label", "centroid", "bbox", "area", 
                  "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length",
                  "min_intensity", "max_intensity", "mean_intensity"]

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
        props_table.insert(props_table.columns.get_loc('label') + 1, column="scale", value=f"{round_scale}")
    else: 
        props_table.insert(props_table.columns.get_loc('label') + 1, column="scale", value=f"{tuple(np.ones(segmentation_img.ndim))}") 

    props_table.insert(props_table.columns.get_loc('volume') + 1, "surface_area", surface_area_tab)
    props_table.insert(props_table.columns.get_loc('surface_area') + 1, "SA_to_volume_ratio", props_table["surface_area"].div(props_table["volume"]))
    props_table[f"{mask_name}_volume"] = mask_vol

    # print this statement to let user known of suppressed warnings
    if Warning: print(f"Warning(s) suppressed while quantifying {seg_name}. See 'method_morphology.ipynb' notebook for more details.")

    return props_table



# quantify the morphology of one or more organelle from one cell
def get_org_morphology(source_file_path: str,
                         list_obj_names: List[str],
                         list_obj_segs: List[np.ndarray],
                         list_intensity_img: List[np.ndarray],
                         list_region_names: List[str],
                         list_region_segs: List[np.ndarray],
                         mask_name: str,
                         scale: Union[tuple,None] = None):
    """
    Measure the amount, size, and shape of multiple organelles from a single cell

    Parameters:
    ----------
    source_file: str
        file path; this is used for recorder keeping of the file name in the output data tables
    list_obj_names: List[str]
        a list of object names (strings) that will be measured; this should match the order in list_obj_segs
    list_obj_segs: List[np.ndarray]
        a list of 3D (ZYX) segmentation np.ndarrays that will be measured per cell; the order should match the list_obj_names 
    list_intensity_img: List[np.ndarray]
        a list of 3D (ZYX) grayscale np.ndarrays that will be used to measure fluoresence intensity in each region and object
    list_region_names: List[str]
        a list of region names (strings); these should include the mask (entire region being measured - usually the cell) 
        and other sub-mask regions from which we can meausure the objects in (ex - nucleus, neurites, soma, etc.). It should 
        also include the centering object used when created the XY distribution bins.
        The order should match the list_region_segs
    list_region_segs: List[np.ndarray]
        a list of 3D (ZYX) binary np.ndarrays of the region masks; the order should match the list_region_names.
    mask: str
        a str of which region name (contained in the list_region_names list) should be used as the main mask (e.g., cell mask)
    scale: Union[tuple,None] = None
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)

    Returns:
    ----------
    Dataframe of measurements of organelle morphology

    """
    print(f"Quantifying organelle morphology from {source_file_path}.")

    # select the mask from the region list
    mask = list_region_segs[list_region_names.index(mask_name)]
    
    # empty list to collect a morphology data for each organelle
    org_tabs = []

    # loop through the list of organelles and run the get_morphology_metrics function
    for j, target in enumerate(list_obj_names):
        # select intensity image
        org_img = list_intensity_img[j]  
        
        # select segmentation and if ER, ensure it is only one object
        if target == 'ER':
            org_obj = (list_obj_segs[j] > 0).astype(np.uint16)  
        else:
            org_obj = list_obj_segs[j]
        
        # run get_morphology_metrics function to output a table of measurements
        org_metrics = get_morphology_metrics(segmentation_img=org_obj, 
                                            seg_name=target,
                                            intensity_img=org_img, 
                                            mask=mask,
                                            mask_name=mask_name,
                                            scale=scale)

        # add table to list above
        org_tabs.append(org_metrics)

    # combine the lists for each organelle into one table
    final_org_tab = pd.concat(org_tabs, ignore_index=True)

    # add a new column to list the name of the image these data are derived from 
    final_org_tab.insert(loc=0,column='image_name',value=source_file_path.stem)
    

    return final_org_tab


# batch process organelle morphology quantification for multiple cells from a single experiment
def batch_process_org_morph(out_file_name: str,
                            seg_path: Union[Path,str],
                            out_path: Union[Path, str], 
                            raw_path: Union[Path,str], 
                            raw_file_type: str,
                            organelle_names: List[str],
                            organelle_channels: List[int],
                            region_names: List[str],
                            mask_name: str,
                            scale:bool=True,
                            seg_suffix:Union[str, None]=None) -> int :
    """  
    batch process segmentation quantification (morphology, distribution, contacts); this function is currently optimized to process images from one file folder per image type (e.g., raw, segmentation)
    the output csv files are saved to the indicated out_path folder

    Parameters:
    ----------
    out_file_name: str
        the prefix to use when naming the output datatables
    seg_path: Union[Path,str]
        Path or str to the folder that contains the segmentation tiff files
    out_path: Union[Path, str]
        Path or str to the folder that the output datatables will be saved to
    raw_path: Union[Path,str]
        Path or str to the folder that contains the raw image files
    raw_file_type: str
        the file type of the raw data; ex - ".tiff", ".czi"
    organelle_names: List[str]
        a list of all organelle names that will be analyzed; the names should be the same as the suffix used to name each of the tiff segmentation files
        Note: the intensity measurements collect per region (from get_region_morphology_3D function) will only be from channels associated to these organelles 
    organelle_channels: List[int]
        a list of channel indices associated to respective organelle staining in the raw image; the indices should listed in same order in which the respective segmentation name is listed in organelle_names
    region_names: List[str]
        a list of regions, or masks, to measure; the order should correlate to the order of the channels in the "masks" output segmentation file
    mask: str
        the name of the region to use as the mask when measuring the organelles; this should be one of the names listed in regions list; usually this will be the "cell" mask
    scale:bool=True
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)
    seg_suffix:Union[str, None]=None
        any additional text that is included in the segmentation tiff files between the file stem and the segmentation suffix

    Returns:
    ----------
    count: int
        the number of images processed
        
    """
    start = time.time()
    count = 0

    # create path objects if inputs are strings
    if isinstance(raw_path, str): raw_path = Path(raw_path)
    if isinstance(seg_path, str): seg_path = Path(seg_path)
    if isinstance(out_path, str): out_path = Path(out_path)
    
    # create directory is it doesn't exist
    if not Path.exists(out_path):
        Path.mkdir(out_path)
        print(f"Output file path not found. Making {out_path}.")
    

    # reading list of files from the raw path
    img_file_list = list_image_files(raw_path, raw_file_type)
    len_file_list = len(img_file_list)

    # list of organelle segmentation and masks files to collect from each image
    segs_to_collect = organelle_names + region_names

    # containers to collect data tabels
    org_tabs = []

    # loop through list of cell analyzing each and appending the data to the empty list
    for img_f in img_file_list:
        img_start = time.time()
        count = count + 1
        filez = find_segmentation_tiff_files(img_f, segs_to_collect, seg_path, seg_suffix)

        # read in raw file and metadata
        img_data, meta_dict = read_czi_image(filez["raw"])

        # create intensities from raw file as list based on the channel order provided
        intensities = [img_data[ch] for ch in organelle_channels]

        # store organelle images as list
        organelles = [read_tiff_image(filez[org]) for org in organelle_names]

        # load regions as a list based on order in list (should match order in "masks" file)
        regions = [read_tiff_image(filez[r]) for r in region_names] 

        # define the scale
        if scale is True:
            scale_tup = meta_dict['scale']
        else:
            scale_tup = None

        org_metrics = get_org_morphology(source_file_path=img_f,
                                            list_obj_names=organelle_names,
                                            list_obj_segs=organelles,
                                            list_intensity_img=intensities, 
                                            list_region_names=region_names,
                                            list_region_segs=regions, 
                                            mask_name=mask_name,
                                            scale=scale_tup)

        org_tabs.append(org_metrics)
        end2 = time.time()
        print(f"Completed quantification of {meta_dict['file_name']} in {(end2-img_start)/60} mins.")
        print(f"{count}/{len_file_list} images have been processed.")
        print(f"Time elapsed: {(end2-start)/60} mins")

    final_org = pd.concat(org_tabs, ignore_index=True)

    org_csv_path = out_path / f"{out_file_name}_org_morph.csv"
    final_org.to_csv(org_csv_path)

    end = time.time()
    print(f"Quantification for {count} files is COMPLETE! Files saved to '{out_path}'.")
    print(f"It took {(end - start)/60} minutes to quantify these files.")
    return final_org


# summarize morphology values per organelle per cell across one or more experiments
def batch_org_morph_summary_stats(csv_path_list: List[str],
                                    out_path: str,
                                    out_preffix: str):
    """" 
    csv_path_list: List[str],
        A list of path strings where .csv files to analyze are located.
    out_path: str,
        A path string where the summary data file will be output to
    out_preffix: str
        The prefix used to name the output file. An "_" will be included between this prefix and the file suffix.
    """
    ds_count = 0
    fl_count = 0
    ###################
    # Read in the csv files and combine them into one of each type
    ###################
    # create empty list to hold the morphology tables from different experiments
    org_tabs = []
    org = "_org_morph"

    # loop through all of the locations listed above and find the _org_morph files; append them to the list above
    for loc in csv_path_list:
        ds_count = ds_count + 1
        files_store = sorted(loc.glob("*.csv"))
        for file in files_store:
            fl_count = fl_count + 1
            stem = file.stem

            if org in stem:
                test_orgs = pd.read_csv(file, index_col=0)
                test_orgs.insert(0, "dataset", stem[:-11])
                org_tabs.append(test_orgs)

    # combine the org_morph lists found above into one table
    org_df = pd.concat(org_tabs,axis=0, join='outer')


    ###################
    # summary stat group
    ###################
    group_by = ['dataset', 'image_name', 'object']
    sharedcolumns = ["SA_to_volume_ratio", "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length"]
    ag_func_standard = ['mean', 'median', 'std']

    ###################
    # summarize shared measurements between org_df and contacts_df
    ###################
    tab1 = org_df[group_by + ['volume']].groupby(group_by).agg(['count', 'sum'] + ag_func_standard)
    tab2 = org_df[group_by + ['surface_area']].groupby(group_by).agg(['sum'] + ag_func_standard)
    tab3 = org_df[group_by + sharedcolumns].groupby(group_by).agg(ag_func_standard)
    org_summary = pd.merge(tab1, tab2, 'outer', on=group_by)
    org_summary = pd.merge(org_summary, tab3, 'outer', on=group_by)
    org_summary[('volume', 'cell')] = org_df[group_by + ['cell_volume']].groupby(group_by).first()


    ###################
    # add normalization
    ###################
    # organelle area fraction
    org_summary[('volume', 'fraction')] = org_summary[('volume', 'sum')]/org_summary[('volume', 'cell')]
    org_summary = org_summary.sort_index(axis=1, level=0, ascending=False)


    ###################
    # flatten datasheets and combine
    ###################
    # org flattening
    org_final = org_summary.unstack(-1)
    new_col_order = ['dataset', 'image_name', 'object', 'volume', 'surface_area', 'SA_to_volume_ratio', 
                    'equivalent_diameter', 'extent', 'euler_number', 'solidity', 'axis_major_length']
    new_cols = org_final.columns.reindex(new_col_order, level=0)
    org_final = org_final.reindex(columns=new_cols[0])
    org_final.columns = ["_".join((col_name[-1], col_name[1], col_name[0])) for col_name in org_final.columns.to_flat_index()]

    #renaming, filling "NaN" with 0 when needed, and removing ER_std columns
    for col in org_final.columns:
        if col.endswith(("_count_volume","_sum_volume", "_mean_volume", "_median_volume")):
            org_final[col] = org_final[col].fillna(0)
        if col.endswith("_count_volume"):
            org_final.rename(columns={col:col.split("_")[0]+"_count"}, inplace=True)
        if col.startswith(("ER_std_", "ER_mean_", "ER_median_")):
            org_final.drop(columns=[col], inplace=True)
    org_final = org_final.reset_index()


    ###################
    # export summary sheets
    ###################
    org_summary.to_csv(str(out_path) + f"/{out_preffix}_per_org_summarystats.csv")


    print(f"Organelle morphology summary for {fl_count} files from {ds_count} dataset(s) is complete.")
    return org_summary