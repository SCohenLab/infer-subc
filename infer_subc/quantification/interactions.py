
import pandas as pd
import numpy as np
import math
from pathlib import Path
from typing import List, Union


from skimage.measure import regionprops_table

from infer_subc.core.img import *
from infer_subc.quantification.stats import *
from infer_subc.quantification.stats_helpers import *
from infer_subc.organelles import * 
from infer_subc.quantification.morphology import get_morphology_metrics
from infer_subc.quantification.csv_io import load_existing_keys_csv, append_atomic_csv
from infer_subc.core.file_io import export_inferred_organelle



def make_dict(list_obj_names: list[str],
              list_obj_segs: list[np.ndarray]) -> dict[str:np.ndarray]:
    '''
    Create a dictionary of organelle segmentations from a list of organelle names and a list of organelle segmentation arrays.

    Parameters
    ----------
    list_obj_names : list[str]
        A list of organelle names as strings.
    list_obj_segs : list[np.ndarray]
        A list of organelle segmentation image arrays as numpy ndarrays.

    Returns
    -------
    organelle_segs : dict[str:np.ndarray]
        A dictionary of organelle segmentations with organelle names as keys and segmentation arrays as values
    '''
    
    organelle_segs = {}
    for idx, name in enumerate(list_obj_names):
        if name == "ER":
            organelle_segs[name]=(list_obj_segs[idx]>0).astype(np.uint8) #if ER, make binary mask (combine all ER objects into one with ID #1)
        else:
            organelle_segs[name]=list_obj_segs[idx]
            
    return organelle_segs


def create_overlap(inter_name:str,
                   organelle_segs: dict[str, np.ndarray],
                   name_splitter: str="X") -> tuple[np.ndarray, np.ndarray]:
    '''
    Create an image of the overlap regions between the selected organelles.

    Parameters
    ----------
    inter_name : str
        A string of organelle names separated by the specified name_splitter.
    organelle_segs : dict[str:np.ndarray]
        A dictionary of organelle segmentations with organelle names as keys and segmentation image arrays as values.
    name_splitter : str, optional
        The character used to split the organelle names in the orgs string, by default "X". 
        For example, "mitoXlyso" would indicate an interaction between mito and lyso.
    
    Returns
    -------
    site : np.ndarray
        An image array of the overlap regions between the selected organelles, with unique integer IDs for each interaction site.
    '''
    site = np.ones_like(organelle_segs[inter_name.split(name_splitter)[0]]) 
    for org in inter_name.split(name_splitter):        
        b = organelle_segs[org]             # select organelle
        valid = (b>0)*(site>0)              # logical and: select MASK of region where site (everything or previous overlap) overlaps with org b
        digit = len(str(np.max(site)))      # find the max ID number in the site image
        site = (b*(10**(digit)))+site       # multiply the organelle b image by 10^digit (up to the next highest multiplier of 10 compared to the original org number per image) then add the site image to that --> preservation of declumped objects
        site[valid.astype(bool)==False]=0   # select everywhere outside the overlap in the site image and get rid of it --> only object multiplied by 10^digit, but now there is separation between declumped neighbors inherited from the original ID numbers added in
        site = label(site)                  # label the overlap to minimize the max number & provide logical number system to contacts in end
    return site


def find_inter_labels(overlap_img: np.ndarray,
                       interaction_name: str,
                       org_dict: dict[str, np.ndarray],
                       name_splitter: str = "X") -> pd.DataFrame:
    '''
    Identify which organelle IDs are involved in each unique interaction site; 
    the organelle ID numbers are joined by underscores and returned in a table of 
    unique identifiers for each interaction site.

    Parameters
    ----------
    overlap_img : np.ndarray
        An image array of the overlap regions between the organelles included in the org_dict variable;
        each interaction site should be labeled with unique integer IDs that will be included in the output table.
    interaction_name : str
        A string of organelle names separated by the specified splitter.
    org_dict : dict[str:np.ndarray]
        A dictionary of organelle segmentations with organelle names as keys and segmentation image arrays as values.
    name_splitter : str, optional
        The character used to split the organelle names in the orgs string, by default "X". 
        For example, "mitoXlyso" would indicate an interaction between mito and lyso.
        
    Returns
    -------
    inter_tab : pd.DataFrame
        A pandas DataFrame table with unique identifiers (integer IDs and labels) associated to each interaction site.
        `ID`: unique integer identifier for each interaction site in the overlap image. Each site will have a different ID number.
        `object`: the name of the interaction sites being examined, created by joining the organelle names with the specified splitter.
        `label`: a string of organelle ID numbers involved in each interaction site, joined by underscores.
    '''

    # use regionprops table to list interaction sites by unique index and extract slice for each object
    props = regionprops_table(overlap_img, properties=['label', 'slice'])

    # create a list of the organelle ID numbers involved in each interaction site
    involved = interaction_name.split(name_splitter)
    indexes = {'ID': [], 'label': []}

    for index, l in enumerate(props["label"]):
        over_inv = []
        for org in involved:
            volume = overlap_img[props["slice"][index]]
            lorg = org_dict[org][props["slice"][index]]
            volume = volume==l
            lorg = lorg[volume]                                 
            all_inv = np.unique(lorg[lorg>0]).tolist()          
            if len(all_inv) != 1:
                print(f"we have an error.  as-> {all_inv}") # ensure that only one org object is 
            over_inv.append(f"{all_inv[0]}")
        indexes['ID'].append(l)
        indexes['label'].append('_'.join(over_inv))

    inter_tab = pd.DataFrame(indexes)
    inter_tab.insert(0, 'object', interaction_name, True)
    
    return inter_tab


def assess_if_higher_order_int(site: np.ndarray,
                                site_name: str,
                                inter_labels_tab: pd.DataFrame,
                                organelle_segs: dict[str:np.ndarray],
                                splitter: str="X"):
    """
    Determine which interaction sites are included in higher order interactions.
    
    An interaction site is considered part of a higher order interaction if it overlaps 
    with an additional organelle not included in the original interaction site definition.
    
    Parameters
    ----------
    site : np.ndarray
        Labeled image of interaction sites with unique integer IDs
    site_name : str
        Interaction site name (organelle names joined by splitter)
    inter_labels_tab : pd.DataFrame
        Table with columns: 'ID', 'object', 'label'
    organelle_segs : dict[str:np.ndarray]
        Dictionary of all organelle segmentations {name: image_array}
    splitter : str, default="X"
        Character separating organelle names in site_name
    
    Returns
    -------
    LOi_NR : np.ndarray
        Image containing only sites NOT in higher-order interactions
    new_tab : pd.DataFrame
        Input table with added 'in_higher_order' boolean column
    """
    
    # Get organelles not in current site
    site_orgs = set(site_name.split(splitter))
    other_orgs = {name: seg for name, seg in organelle_segs.items() if name not in site_orgs}
    
    if not other_orgs:
        # No other organelles to check, all sites are non-redundant
        in_higher_order = pd.Series([False] * len(inter_labels_tab))
        new_tab = inter_labels_tab.copy()
        new_tab.insert(inter_labels_tab.columns.get_loc('label') + 1, "in_higher_order", in_higher_order)
        return site.copy(), new_tab
    
    # Create combined mask of all other organelles
    other_orgs_mask = np.zeros_like(site, dtype=bool)
    for org_seg in other_orgs.values():
        other_orgs_mask |= (org_seg > 0)
    
    # Find sites that overlap with any other organelle
    overlap_mask = (site > 0) & other_orgs_mask
    higher_order_ids = set(np.unique(site[overlap_mask])) - {0}
    
    # Create output image (zero out higher-order sites)
    LOi_NR = site.copy()
    if higher_order_ids:
        mask_to_remove = np.isin(site, list(higher_order_ids))
        LOi_NR[mask_to_remove] = 0
    
    # Mark higher-order sites in table
    in_higher_order = inter_labels_tab['ID'].isin(higher_order_ids)
    new_tab = inter_labels_tab.copy()
    new_tab.insert(inter_labels_tab.columns.get_loc('label') + 1, "in_higher_order", in_higher_order)
    
    return LOi_NR, new_tab


def create_interaction_sites(interaction_orgs: List[str], 
                              org_name_list:List[str],
                              org_seg_list: List[np.ndarray],
                              name_splitter: str="X",
                              mask: Union[np.ndarray, None]=None,
                              mask_name: Union[str, None]=None) -> tuple[np.ndarray, np.ndarray, pd.DataFrame]:
    
    '''
    Create an image of the overlap regions between the selected organelles and a table of unique identifiers 
    associated to each interaction site. The entire image is used to create the interaction sites; no mask is applied.
    
    Parameters 
    ----------
    interaction_orgs : List[str]
        A list of organelle names as strings to be included in the interaction site.
    org_name_list : List[str]
        A list of all organelle names as strings. These will be the organelles used to create interaction sites.
    org_seg_list : List[np.ndarray]
        A list of all organelle segmentation images as numpy ndarrays. These should be in the same order as the 
        org_name_list list.
    name_splitter : str, optional
        The character used to separate the organelle names in the org_name_list string within the new interaction 
        site name, by default "X". For example, "mitoXlyso" would indicate an interaction between mito and lyso. 
        Use of other splitters may cause issues during downstream analysis in infer-subc. 
        Specifically, avoid using "_" or "-" as a splitter as they are used in other parts of the analysis.
    mask : Union[np.ndarray, None], optional
        A binary np.ndarray mask of the area to measure from. If None, the whole image is analyzed.
    mask_name : Union[str, None], optional
        The name of the mask region being analyzed. If None, the whole image is analyzed.

        
    Returns
    -------
    overlap_img : np.ndarray
        An image array of the overlap regions between the selected organelles, with unique integer IDs for each 
        interaction site.   
    lower_order_sites : np.ndarray
        An image array of the interaction sites that are NOT involved in higher order interactions
    inter_tab : pd.DataFrame
        A pandas DataFrame table with unique identifiers (integer IDs and labels) associated to each interaction site.
        `ID`: unique integer identifier for each interaction site in the overlap image. Each site will have a 
        different ID number.
        `object`: the name of the interaction sites being examined, created by joining the organelle names with the 
        specified splitter.
        `label`: a string of organelle ID numbers involved in each interaction site, joined by underscores.
        TODO: consider changing label "splitter" to "X" to be consistent with the name_splitter.
    '''
    
    if len(interaction_orgs)<2:
        raise ValueError("Please select at least two organelles to define an interaction site.")
    else:
        # create name for interaction site of just organelles to include in the interaction site
        interaction_name = name_splitter.join(interaction_orgs)

        # run function to create dictionary of organelle segmentations - including all organelles
        org_dict = make_dict(org_name_list, org_seg_list)

        # create the overlap image
        overlap_img = create_overlap(interaction_name, org_dict)

        # apply mask to overlap image if provided
        if mask_name is None and mask is not None:
            raise ValueError("The mask_name parameter must be provided if mask is not None")
        elif mask is None and mask_name is None:
            input_labels = overlap_img
            mask_name = "whole_image"
        else:
            input_labels = label(apply_mask(overlap_img, mask)).astype(int)

        # use regionprops table to list interaction sites by unique index and extract slice for each object
        inter_tab = find_inter_labels(input_labels, interaction_name, org_dict)
        inter_tab.insert(0, 'mask_name', mask_name, True)

        # determine if each site is also involved in a higher order interaction (there are more than the specified organelles involved)
        lower_order_sites, inter_tab = assess_if_higher_order_int(input_labels, interaction_name, inter_tab, org_dict)

        return input_labels, lower_order_sites, inter_tab
    

def create_interaction_degrees(org_name_list:List[str],
                               org_seg_list: List[np.ndarray],
                               mask: Union[np.ndarray, None]=None,
                               mask_name: Union[str, None]=None,
                               scale: Union[tuple, None]=None) -> (np.ndarray, pd.DataFrame):
    '''
    Create interaction degree image and quantification table for a set of organelle segmentations within a cell mask.

    Parameters
    ----------
    org_name_list : List[str]
        A list of organelle names as strings. These will be the organelles used to create interaction sites.
    org_seg_list : List[np.ndarray]
        A list of organelle segmentation images as numpy ndarrays. These should be in the same order as the org_name_list list.
    mask : np.ndarray
        A binary image array representing the cell (or other) mask. 
        If None is provided, a whole image mask will be used.
    mask_name : str, optional
        The name of the cell mask, by default "cell".
        If None, "whole_image" will be used.
    scale : Union[tuple, None], optional
        A tuple representing the scale of the image in ZYX dimensions, by default None. 
        If None, a scale of (1,1,1) will be used.   
    
    Returns
    -------
    all_orgs : np.ndarray
        An image array representing the degree of interactions between the organelles.
        The cell mask was not applied to this image, but was applied prior to quantification below.
        Each voxel value indicates the number of organelles present at that location.
    final_quant_tab : pd.DataFrame
        A pandas DataFrame table with quantification of the interaction degrees for each organelle and the cell mask.
        The table includes voxel counts and volumes for each degree of interaction.
    '''
    # add all binary organelle segmentation masks together into the new all_orgs object
    all_orgs = np.zeros_like(org_seg_list[0], dtype=np.uint8)
    for o in org_seg_list:
        all_orgs = all_orgs + (o>0)
    
    # create empty dictionary to hold quantification results
    quant_tabs = [] 

    # fill in a scale value if none is specified
    if scale is None:
        scale = (1,1,1) 

    # safe guard against no mask being provided
    if mask is None:
        mask = np.ones_like(org_seg_list[0], dtype=bool)
        mask_name = "whole image"
    
    # loop through each organelle and the cell mask
    for reg_name, reg in zip(([mask_name]+org_name_list), ([mask]+org_seg_list)):
        # mask with cell mask and then with region of interest
        masked = apply_mask(reg, mask.astype(bool))
        all_orgs_masked = apply_mask(all_orgs, masked.astype(bool))

        # count the number of voxels with each degree of interaction per cell
        degrees, counts = np.unique(all_orgs_masked, return_counts=True)
        nway_quant = dict(zip(degrees, counts))
        nway_quant.pop(0, None)  # remove background count

        # calculate additional metrics
        tot_org_vox = np.sum(list(nway_quant.values()))
        tot_reg_vox = np.count_nonzero(masked>0)

        # create dictionary of unscaled results
        nway_quant = {'voxel_count_with_0_org(s)': tot_reg_vox - tot_org_vox,
                      **{"voxel_count_with_" + str(key) + "_orgs(s)": value for key, value in nway_quant.items()},
                      'voxel_count_region': tot_reg_vox}

        # created scaled dictionary
        prod_scale = math.prod(scale)
        nway_quant_scaled = {k.replace('voxel_count', 'volume'): v * prod_scale for k, v in nway_quant.items() if 'voxel_count' in k}

        # combinde unscaled and scale results
        final_dict = {'scale': (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4)),
                      'object': reg_name,
                      **nway_quant, **nway_quant_scaled}

        # make it into a dataframe
        quant_tabs.append(final_dict)

    # combine into one table
    final_quant_tab = pd.DataFrame(quant_tabs)
    final_quant_tab.insert(0, column="mask_name", value=mask_name)

    return all_orgs, final_quant_tab


def all_combos(list_obj_names: list[str], 
               splitter: str="X") -> list:
    """
    Create names for all possible combinations of organelle interaction site types from a list of organelles

    Parameters
    ----------
    list_obj_names: list[str], 
        a list of names as strings for the organelles segmented in the image being analyzed; organelle names should match the naming suffix on the organelle segmentation file
        ex) mitochondria file name: "img1-mito.tiff" 
            naming suffix: "mito"
            list of organelles: ["mito", "lyso", "perox", ...]
    splitter: str="X"
        a character you wish to use as the seperator between organelle names when creating interaction site names
        "X" is the recommended splitter

    Output
    ------
    possib: dict
        a list of the names for all possible organelle interaction site combinations

    """
    all_pos = []
    for n in list(map(lambda x:x+2, (range(len(list_obj_names)-1)))):
        all_pos += itertools.combinations(list_obj_names, n)
    possib = [splitter.join(inter) for inter in all_pos]
    return possib


def get_interaction_metrics(source_file_path: str,
                             list_obj_names: List[str],
                             list_obj_segs: List[np.ndarray],
                             list_intensity_img: Union[List[np.ndarray], None]=None,
                             list_region_names: Union[List[str], None]=None,
                             list_region_segs: Union[List[np.ndarray], None]=None,
                             mask_name: Union[str, None]=None,
                             scale: Union[tuple, None]=None,
                             splitter: str="X",
                             include_morpho:bool=True,
                             include_inter_degrees:bool=True,
                             include_dist:bool=True, 
                             dist_centering_obj: Union[str, None]=None,
                             dist_num_bins: Union[int, None]=5,
                             dist_center_on: Union[bool, None]=False,
                             dist_keep_center_as_bin: Union[bool, None]=True,
                             dist_zernike_degrees: Union[int, None]=9) -> tuple[dict, pd.DataFrame, Union[pd.DataFrame, None], 
                                                                                Union[np.ndarray, None], Union[pd.DataFrame, None], 
                                                                                Union[np.ndarray, None], Union[np.ndarray, None]]:
   
    """
    Quantify organelle interaction metrics including morphology, distribution, and degree of interactions for a image or region
    (e.g., cell) within an image.
    
    Parameters
    ----------
    source_file_path : str or Path
        Path to the source image file. This will be used as part of the metadata information in the output tables. 
        The input images are not derived from this path, but rather are provided directly as arrays in the list_obj_segs and 
        list_intensity_img variables below.
    list_obj_names : List[str]
        List of organelle names. These names should match the suffix on the segmentation image files.
    list_obj_segs : List[np.ndarray]
        List of 3D organelle segmentation arrays matching the order included in list_obj_names.
    list_intensity_img : Union[List[np.ndarray], None]=None
        List of 3D intensity channels from the raw image used to produce the segmentations in list_obj_segs.
        The order here should match the list_obj_segs and list_obj_names variables.
        Additional intensity channels not matching one of the segmented organelles/included in list_obj_names should not be included.
        If no intensity analysis is to be included, specify None here.
    list_region_names : Union[List[str], None]=None
        List of segmented region/mask names. These names should match the suffix on the segmentation image files.
        This should include:
            - a mask segmentation, such as the cell mask, for masking during all interactions analysis; else, the entire image will be 
            quantified. Only one objects per mask image will be analyzed. If there are more than one included, they will be combined 
            prior to analysis and the entire region will be quantified. If no mask is provided, the entire image will be quantified.
            - a centering object, such as the nucleus, for distribution analysis; else the center of the mask region will be used as 
            the XY distribution centering point if distribution analysis is included.
    list_region_segs : Union[List[np.ndarray], None]=None
        List of 3D region segmentation arrays matching the order specified in list_region_names. Specify None if no regions are provided.
    mask_name : Union[str, None]=None
        Name of the region to use as the mask for analysis; if not specified, the entire image will be quantified.
    splitter : str, default="X"
        Character used to separate organelles within the interaction site names
        Ex) "mitoXlyso" for mito-lyso interactions
    scale : Union[tuple, None], default=None
        Name of the region to use as the mask for analysis; if not specified, the entire image will be quantified.
    include_morpho : bool, default=True
        Whether to compute morphology metrics for each interaction site.
    channel_axis : int, default=0
        The index of the channel dimension axis in the intensity image.
    include_inter_degrees : bool, default=True
        Whether to compute interaction degree analysis for the entire image or mask region.
    include_dist : bool, default=True
        Whether to compute distribution metrics for the each interaction site type.
    dist_centering_obj : Union[str, None], default=None
        Name of the region to use for centering distribution analysis.
        This region should be included in the list_region_names and list_region_segs variables.
        If not specified, the center of the mask, or entire image if no mask was specified, will be used as the centering object.
    dist_num_bins : Union[int, None], default=5
        Number of radial bins to create in the XY distribution analysis.
    dist_center_on : Union[bool, None], default=True
        Whether to start creation of the XY bins from the center (True) or the edge (False) of the centering object.
    dist_keep_center_as_bin : Union[bool, None], default=True
        Whether to keep the centering object as the first XY bin. 
    dist_zernike_degrees : Union[int, None], default=9
        Zernike polynomial degree for circular shape/pattern analysis in the XY distribution analysis.
        If None and include_dist=True, no Zernike features will be calculated.
    
    Returns
    -------
    inter_sites : dict
        Dictionary of interaction site name, np.ndarray image pairs for all possible interaction site combinations in the cell
    morph_final_combo : pd.DataFrame
        Combined morphology metrics for all interaction sites of each interaction type.
        If include_morpho=False, this will only include the interaction site metrics calculated in the
        infer_subc.quantification.interactions.create_interaction_sites() function, including which organelles are involved in each site
        and if they are in higher order interaction sites; it will not list morphology metrics for each interaction site.
    dist_final_combo : pd.DataFrame or None
        XY and Z distribution metrics for each interaction site (if include_dist=True)
    degree_img : np.ndarray or None
        Degree of interactions image (if include_inter_degrees=True)
    degree_tab : pd.DataFrame or None
        Degree of interactions table (if include_inter_degrees=True)
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
    print(f"Quantifying organelle interactions from {source_file_path.name}")

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
        mask = (list_region_segs[list_region_names.index(mask_name)] > 0).astype(int) # ensure mask is binary and integer type for later multiplication with segmentation images
        print(f"Mask '{mask_name}' will be applied before analysis.")

    # list all possible interaction site types based on the org_file_names list specified above
    possib_int_types = all_combos(list_obj_names, splitter)

    # recreate raw_intensity image based on list intensity channels above to ensure proper order
    if include_morpho:
        if list_intensity_img is None:
            intensity_img = None
            print("No intensity images provided. Morphology metrics that require intensity images will not be calculated.")
        else:
            intensity_img = np.stack(list_intensity_img)

    # collect centering object image
    if include_dist:
        if dist_centering_obj is None:
            print("No centering object provided. Using center of mask or entire image for distribution centering.")
            centering_img = None
        elif dist_centering_obj not in list_region_names:
            raise ValueError(f"Centering object '{dist_centering_obj}' not found in region names: {list_region_names}")
        else:
            centering_img = list_region_segs[list_region_names.index(dist_centering_obj)]


    # collect interaction metric tabs
    morph_combo_tabs = []
    dist_combo_tabs = []
    XY_bins_imgs = []
    XY_wedges_imgs = []
    inter_sites = {}

    # loop through interaction site types and create interaction sites, measure morphology and distributions
    for overlap_ID in possib_int_types:
        # list organelles included in this interaction site only
        orgs_included = overlap_ID.split("X")

        # create interaction site & metadata information
        inter_obj, lower_ord_sites, inter_tab = create_interaction_sites(orgs_included,
                                                                         list_obj_names,
                                                                         list_obj_segs, 
                                                                         name_splitter=splitter,
                                                                         mask=mask,
                                                                         mask_name=mask_name)
        del lower_ord_sites
        inter_sites[overlap_ID] = inter_obj

        # measure interaction site morphology
        if include_morpho:
            morpho_metrics = get_morphology_metrics(segmentation_img=inter_obj, 
                                                    seg_name=overlap_ID,
                                                    intensity_img=intensity_img, 
                                                    intensity_ch_names=list_obj_names,
                                                    channel_axis=0,
                                                    mask=mask,
                                                    mask_name=mask_name,
                                                    scale=scale)
            morpho_metrics.rename(columns={'label':'ID'}, inplace=True)
            inter_tab = pd.merge(inter_tab, morpho_metrics, how='right', on=['object', 'ID', 'mask_name'])
            inter_tab['in_higher_order'] = inter_tab['in_higher_order'].astype(bool) # force to boolean type
        
        morph_combo_tabs.append(inter_tab)

        # measure interaction site distibutions
        if include_dist:
            XY_distribution, XY_bins, XY_wedges = get_XY_distribution(mask=mask,
                                                                      mask_name=mask_name,
                                                                    centering_obj=centering_img,
                                                                    obj=inter_obj,
                                                                    obj_name=overlap_ID,
                                                                    scale=scale,
                                                                    num_bins=dist_num_bins,
                                                                    center_on=dist_center_on,
                                                                    keep_center_as_bin=dist_keep_center_as_bin,
                                                                    zernike_degrees=dist_zernike_degrees)
            # if XY_bins_imgs list is empty append, if not skip
            if not XY_bins_imgs and not XY_wedges_imgs:
                XY_bins_imgs.append(XY_bins)
                XY_wedges_imgs.append(XY_wedges)

            Z_distribution = get_Z_distribution(mask=mask, 
                                                mask_name=mask_name,
                                                obj=inter_obj,
                                                obj_name=overlap_ID,
                                                center_obj=centering_img,
                                                scale=scale)

            interaction_dist_tab = pd.merge(XY_distribution, Z_distribution)
            dist_combo_tabs.append(interaction_dist_tab)
            
    # merge the tables together
    morph_final_combo = pd.concat(morph_combo_tabs)
    morph_final_combo.insert(loc=0,column='image_name',value=source_file_path.stem)

    if include_dist:
        dist_final_combo = pd.concat(dist_combo_tabs)
        dist_final_combo.insert(loc=0,column='image_name',value=source_file_path.stem)
    else:
        dist_final_combo = None
        XY_bins_imgs = [None]
        XY_wedges_imgs = [None]

    if include_inter_degrees:
        degree_img, degree_tab = create_interaction_degrees(org_name_list=list_obj_names,
                                                            org_seg_list=list_obj_segs,
                                                            mask=mask,
                                                            mask_name=mask_name,
                                                            scale=scale)

        # add source image name to degree table
        degree_tab.insert(loc=0,column='image_name',value=source_file_path.stem)
    else:
        degree_img = None
        degree_tab = None
        
    return inter_sites, morph_final_combo, dist_final_combo, degree_img, degree_tab, XY_bins_imgs[0], XY_wedges_imgs[0] 



def batch_process_interactions_quant(dataset_name: str,
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
                                      int_splitter:str="X",
                                      include_morpho:bool=True,
                                      include_inter_degrees:bool=True,
                                      include_dist:bool=True, 
                                      dist_centering_obj: Union[str, None]=None,
                                      dist_num_bins: Union[int, None]=5,
                                      dist_center_on: Union[bool, None]=False,
                                      dist_keep_center_as_bin: Union[bool, None]=True,
                                      dist_zernike_degrees: Union[int, None]=9,
                                      export_inter_degree_imgs:bool=True,
                                      export_interaction_sites:bool=True,
                                      export_distribution_bins_imgs:bool=True) -> None:
    """
    Batch process interaction quantification for a single dataset (e.g., images collected on the same data). 
    Morphology, distribution, and degree of interaction metrics analysis are all optionally available. 
    Interaction site segmentations and degree of interaction images can also be exported.
    
    Parameters:
    -----------
    dataset_name : str
        A unique string identifier for the dataset being processed. It will be included as metadata in output tables and as 
        part of the output files names. It will be used to identify if any data has already been collected for this dataset.
    raw_path : Union[Path,str]
        Path or str to the folder that contains the raw image files.
    seg_path : Union[Path,str]
        Path or str to the folder that contains the segmentation tiff files.
    quant_path : Union[Path, str]
        Path or str to the folder that the output datatables will be saved to.
    raw_file_type : str
        File type of the raw images (e.g., "czi", "tiff")
    channel_axis : int
        Axis corresponding to the channels in the image data
    organelle_names : List[str]
        List of organelle names to analyze. These names should match the suffix on the organelle segmentation files
    organelle_channels : Union[List[int], None]
        List of intensity channel indices in the raw files corresponding to each organelle included in organelle_names.
        The order should match organelle_names. 
        If no intensity analysis is to be included, specify None here.
    region_names : Union[List[str], None]
        List of region names to analyze. Usually ['cell', 'nuc'] for cell mask and nucleus.
        If no regions are to be included, specify None here.
    mask_name : Union[str, None]
        Name of the mask to use for segmentation (if any). This name should be included in the regions_name variable.
        If None, the entire image will be quantified.
    use_scale : bool
        Whether to apply scaling to the quantitative data; scaled data will be in real world units (e.g., microns) rather than pixels/voxels
    seg_suffix : Union[str, None]
        Any additional text that is included in the segmentation tiff files between the file stem and the segmentation suffix, not including the initial "-"
    int_splitter: str
        Character used to separate organelles within the interaction site names
        Ex) "mitoXlyso" for mito-lyso interactions
        include_morpho : bool, default=True
        Whether to compute morphology metrics
    include_morpho : bool, default=True
        Whether to compute morphology metrics for each interaction site
    include_inter_degrees : bool, default=True
        Whether to compute interaction degree analysis
    include_dist : bool, default=True
        Whether to compute distribution metrics
    dist_centering_obj : str or None, default=None
        Name of the region to use for centering distribution analysis
        This region should be included in the list_region_names and list_region_segs variables
        If not specified, the center of the mask, or entire image if no mask was specified, will be used as the centering object
    dist_num_bins : int or None, default=5
        Number of radial bins to create in the XY distribution analysis
    dist_center_on : bool or None, default=True
        Whether to start creation of the XY bins from the center of the centering object (True) or edge (False)
    dist_keep_center_as_bin : bool or None
        Whether to keep centering object as the first XY bin
    dist_zernike_degrees : int or None, default=9
        Zernike polynomial degree for shape analysis in the XY distribution analysis
        If None and include_dist=True, no Zernike features will be calculated
    export_inter_degree_imgs : bool
        Whether to export interaction degree images
    export_interaction_sites : bool
        Whether to export interaction site images (including interaction site objects across the entire image; not masked)
    export_distribution_bins_imgs : bool
        Whether to export the XY distribution bins and wedges images

    Returns:
    --------
    None
        Saves output files to the specified quantification path
    """

    batch_start = time.time()
    count = 0

    # format/make file paths
    if isinstance(raw_path, str): raw_path = Path(raw_path)
    if isinstance(seg_path, str): seg_path = Path(seg_path)
    if isinstance(quant_path, str): quant_path = Path(quant_path)

    if not Path.exists(quant_path):
        Path.mkdir(quant_path)
        print(f"making {quant_path}")


    # check if any existing data is present in outfiles to skip already processed images
    unique_keys = ['dataset', 'image_name']

    if include_morpho:
        morpho_tab_path = quant_path / f"{dataset_name}-interactions_morphology_metrics.csv"
        existing_morpho_keys = load_existing_keys_csv(morpho_tab_path, unique_keys)
    elif not include_morpho:
        morpho_tab_path = quant_path / f"{dataset_name}-interactions_labels.csv"
        existing_morpho_keys = load_existing_keys_csv(morpho_tab_path, unique_keys)

    if include_dist:
        dist_tab_path = quant_path / f"{dataset_name}-interactions_distribution_metrics.csv"
        existing_dist_keys = load_existing_keys_csv(dist_tab_path, unique_keys)
    else:
        existing_dist_keys = set()

    if include_inter_degrees:
        int_degree_tab_path = quant_path / f"{dataset_name}-interactions_degree_metrics.csv"
        existing_int_degree_keys = load_existing_keys_csv(int_degree_tab_path, unique_keys)
    else:
        existing_int_degree_keys = set()

    if existing_morpho_keys == existing_dist_keys == existing_int_degree_keys:
        existing_keys = existing_morpho_keys
    else:
        existing_keys = existing_morpho_keys.intersection(existing_dist_keys).intersection(existing_int_degree_keys)

    int_degree_img_path = quant_path / f"{dataset_name}-interaction_degree_images" if export_inter_degree_imgs else None
    interaction_sites_path = quant_path / f"{dataset_name}-interaction_site_segmentations" if export_interaction_sites else None
    dist_bins_path = quant_path / f"{dataset_name}-distribution_bins_imgs" if export_distribution_bins_imgs else None

    # list of organelle segmentation and masks files to collect from each image
    segs_to_collect = organelle_names + region_names if region_names is not None else organelle_names

    # reading list of files from the raw path
    img_file_list = list_image_files(raw_path, raw_file_type)
    len_file_list = len(img_file_list)

    # loop through list of cell analyzing each and appending the data to the empty list
    for img_f in img_file_list:
        img_start = time.time()
        count = count+1
        # skip files that have already been processed
        if (dataset_name, img_f.stem) in existing_keys:
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

            # load regions as a list based on order in list
            if region_names is None:
                regions = None
            else:
                regions = [read_tiff_image(filez[r]) for r in region_names]

            # define the scale
            if use_scale is True:
                scale_tup = meta_dict['scale']
            else:
                scale_tup = None

            inter_dict, morph_tab, dist_tab, int_degree_img, int_degree_tab, XY_bins, XY_wedges = get_interaction_metrics(source_file_path=img_f,
                                                                                                    list_obj_names=organelle_names,
                                                                                                    list_obj_segs=organelles,
                                                                                                    list_intensity_img=intensities,
                                                                                                    list_region_names=region_names,
                                                                                                    list_region_segs=regions,
                                                                                                    mask_name=mask_name,
                                                                                                    splitter=int_splitter,
                                                                                                    scale=scale_tup,
                                                                                                    include_morpho=include_morpho,
                                                                                                    include_inter_degrees=include_inter_degrees,
                                                                                                    include_dist=include_dist, 
                                                                                                    dist_centering_obj=dist_centering_obj,
                                                                                                    dist_num_bins=dist_num_bins,
                                                                                                    dist_center_on=dist_center_on,
                                                                                                    dist_keep_center_as_bin=dist_keep_center_as_bin,
                                                                                                    dist_zernike_degrees=dist_zernike_degrees)

            # save the morphology table data per image directly to csv
            if (dataset_name, img_f.stem) in existing_morpho_keys:
                print(f"Skipping morphology metrics for {img_f.name} as it is already listed in the output file.")
            else:
                morph_tab.insert(loc=0,column='dataset',value=dataset_name)
                append_atomic_csv(morpho_tab_path, morph_tab)
            del morph_tab  # free up memory

            # save the interaction site images
            if export_interaction_sites:
                inter_site_cnt=0
                for inter_name, inter_img in inter_dict.items():
                    inter_site_cnt+=1
                    if not (Path(interaction_sites_path)/f"{meta_dict['file_name'].stem}-{inter_name}.tiff").exists():
                        export_inferred_organelle(inter_img.astype(np.uint16), f"{inter_name}", meta_dict, interaction_sites_path)
                    else:
                        if inter_site_cnt<=1:
                            warnings.warn(f"Some of the interaction site images already exist for {meta_dict['file_name'].stem} in {interaction_sites_path}. They will not be overwritten.", UserWarning)
            del inter_dict  # free up memory
           
            # save the distribution table data per image directly to csv
            if include_dist:
                if (dataset_name, img_f.stem) in existing_dist_keys:
                    print(f"Skipping distribution metrics for {img_f.name} as it is already listed in the output file.")
                else:
                    dist_tab = dist_tab.astype(str)  # ensure all data is string to avoid dtype issues
                    dist_tab.insert(loc=0,column='dataset',value=dataset_name)
                    append_atomic_csv(dist_tab_path, dist_tab)

                    if export_distribution_bins_imgs:
                        # export XY bins and wedges as images
                        if not Path(dist_bins_path / f"{img_f.stem}-XY_bins.tiff").exists():
                            export_inferred_organelle(XY_bins.astype(np.uint16), "XY_bins", meta_dict, dist_bins_path)
                        else:
                            warnings.warn(f"The XY distribution bins images already exist for {meta_dict['file_name'].stem} in {dist_bins_path}. They will not be overwritten.", UserWarning)

                        if not Path(dist_bins_path / f"{img_f.stem}-XY_wedges.tiff").exists():
                            export_inferred_organelle(XY_wedges.astype(np.uint16), "XY_wedges", meta_dict, dist_bins_path)
                        else:
                            warnings.warn(f"The XY distribution wedges images already exist for {meta_dict['file_name'].stem} in {dist_bins_path}. They will not be overwritten.", UserWarning)
            del dist_tab  # free up memory
            del XY_bins  # free up memory
            del XY_wedges  # free up memory

            # save the degree table data per image directly to csv
            if include_inter_degrees:
                if (dataset_name, img_f.stem) in existing_int_degree_keys:
                    print(f"Skipping interaction degree metrics for {img_f.name} as it is already listed in the output file.")
                else:
                    int_degree_tab.insert(loc=0,column='dataset',value=dataset_name)
                    append_atomic_csv(int_degree_tab_path, int_degree_tab)

                # save the degree image
                if export_inter_degree_imgs:
                    if not (Path(int_degree_img_path)/f"{img_f.stem}-interactions_degree.tiff").exists():
                        export_inferred_organelle(int_degree_img.astype(np.uint16), f"interactions_degree", meta_dict, int_degree_img_path)
                    else:
                        warnings.warn(f"The {img_f.stem}-interactions_degree.tiff image already exists in {int_degree_img_path}. It will not be overwritten.")
            del int_degree_tab  # free up memory
            del int_degree_img  # free up memory
            
            end2 = time.time()
            print(f"Completed quantification of {img_f.name} in {(end2-img_start)/60} mins.")
            print(f"{count}/{len_file_list} images have been processed.")
            print(f"Time elapsed: {(end2-img_start)/60} mins")

    batch_end = time.time()
    print(f"Quantification for {count} files is COMPLETE! Files saved to '{quant_path}'.")
    print(f"It took {(batch_end - batch_start)/60} minutes to quantify these files.")



def perorg_interactions_cnt(interaction_morpho_df:pd.DataFrame, 
                             org_list:List[str],
                             splitter:str="X") -> pd.DataFrame:
    """
    Summarize interaction counts and volumes per organelle object from interaction morphology data.
    
    Transforms interaction site data (e.g., "mitoXER" with label "06_01") into per-organelle 
    summaries showing how many times each organelle participates in different interaction types.
    
    Parameters
    ----------
    interaction_morpho_df : pd.DataFrame
        The interactions morphology dataframe created by infer_subc.quantification.interactions.batch_process_interactions_quant() 
        or infer_subc.quantification.interactions.get_interaction_metrics() functions.
        The dataFrame must containing the following columns:
        - dataset: experiment identifier
        - image_name: cell/image identifier
        - mask_name: mask identifier
        - scale: image scale information
        - object: interaction site name (e.g., "mitoXER", "mitoXlysoXgolgi")
        - ID: interaction site ID
        - label: underscore-separated organelle IDs (e.g., "06_01")
        - volume: interaction site volume
    org_list : List[str]
        List of all organelle names included in the interactions analysis
    splitter : str, default="X"
        Character used to split interaction site names
    
    Returns
    -------
    pd.DataFrame
        Per-organelle summary with columns:
        - dataset, image_name, object, label
        - num_interaction_types: the number of different interaction types per organelle objects
        - {interaction_type}_count: number of sites of each interaction type (frequency of each interaction type) per organelle objects
        - {interaction_type}_volume: total volume of each interaction type per organelle objects
    """

    # Select and copy data
    meta_cols = ["dataset", "image_name", "mask_name", "scale", "object", "label"]
    df = interaction_morpho_df[meta_cols + ["volume"]].copy()
    
    # # Split columns
    df[['orgs', 'ids']] = df.apply(lambda row: pd.Series([row['object'].split(splitter), row['label'].split('_')]),axis=1)
    
    # Explode to create one row per organelle in each interaction
    records = []
    for _, row in df.iterrows():
        for org, org_id in zip(row['orgs'], row['ids']):
            records.append({'dataset': row['dataset'],
                            'image_name': row['image_name'],
                            'mask_name': row['mask_name'],
                            'scale': row['scale'],
                            'object': org,
                            'label': int(org_id),
                            'interaction_type': row['object'],
                            'volume': row['volume']})

    expanded = pd.DataFrame(records)

    # Summarize interaction sites per organelle object
    agg_dict = {'volume': ['count', 'sum']}
    grouped = expanded.groupby(meta_cols + ['interaction_type']).agg(agg_dict)
    
    grouped.columns = ['count', 'volume']
    grouped = grouped.reset_index()

    # Add interaction degree
    num_inter_types = grouped.groupby(meta_cols)['interaction_type'].nunique().reset_index(name='num_interaction_types')
    
    # Pivot to wide format
    count_pivot = grouped.pivot_table(index=meta_cols, 
                                      columns='interaction_type',
                                      values='count',
                                      fill_value=0).add_suffix('_count')
    
    volume_pivot = grouped.pivot_table(index=meta_cols,
                                       columns='interaction_type',
                                       values='volume',
                                       fill_value=0).add_suffix('_volume')
    
    # Combine
    result = pd.concat([count_pivot, volume_pivot], axis=1).reset_index()
    combo = pd.merge(num_inter_types, result, on=meta_cols)
    
    # Ensure all interaction types present
    all_possible = all_combos(org_list, splitter=splitter)
    for interaction_type in all_possible:
        if f"{interaction_type}_count" not in combo.columns:
            combo[f"{interaction_type}_count"] = 0
        if f"{interaction_type}_volume" not in combo.columns:
            combo[f"{interaction_type}_volume"] = 0
    
    combo['label'] = combo['label'].astype("Int64")

    # fill NA with 0 and format to float values
    num_cols = [col for col in list(combo.columns) if col not in set(meta_cols)]
    combo[num_cols] = combo[num_cols].fillna(0).astype(float)
    
    return combo


def batch_interactions_summary_stats(out_prefix: str,
                                      csv_path_list: List[str],
                                      out_path: str,
                                      organelle_names: List[str],
                                      mask_name: Union[str, None],
                                      splitter: str = "X"):
    """" 
    Batch process interaction quantification summary statistics from multiple datasets.

    Parameters:
    -----------
    out_prefix: str
        The prefix used to name the output file. An "_" will be included between this prefix and the file suffix.
    csv_path_list: List[str],
        A list of path strings where .csv files to analyze are located.
    out_path: str,
        A path string where the summary data file will be output to
    organelle_names: List[str],
        A list of organelle names used in the interaction quantification analysis.
    mask_name: Union[str, None],
        The name of the mask to be used for filtering interaction data.
    splitter: str, default="X"
        The character used to split interaction site names.
    """
    # for keeping track of dataset and file numbers
    ds_count = 0
    fl_count = 0

    ###############################################################
    # Read in the csv files and combine them into one of each type
    ###############################################################
    # create empty list to hold the morphology tables from different experiments
    int_labs = []
    int_morph = []
    int_dist = []
    int_degree = []

    # loop through all of the locations listed above and find the _org_morph files; append them to the list above
    for loc in csv_path_list:
        # list all csv files in the location
        files_store = sorted(loc.glob("*.csv"))

        # find the unique datasets in this location based on the prefixes before "-interactions_"
        prefixes = set(f.name.split("-interactions_")[0] for f in files_store)
        print(f"Found the following datasets in {loc}:", prefixes)
        for prefix in prefixes:
            ds_count = ds_count + 1
            files_subset = [f for f in files_store if f.name.startswith(prefix +"-interactions")]

            # if both morphology and labels files are present, remove the labels file from the list to be processed
            if any("-interactions_morphology_metrics.csv" in f.name for f in files_subset) and any("-interactions_labels.csv" in f.name for f in files_subset):
                    files_subset = [f for f in files_subset if not "-interactions_labels.csv" in f.name]

            for file in files_subset:
                fl_count = fl_count + 1
                stem = file.stem
                
                if "-interactions_labels" in stem:
                    inter_labels = pd.read_csv(file)
                    int_labs.append(inter_labels)
                elif "-interactions_morphology_metrics" in stem:
                    morph = pd.read_csv(file)
                    int_morph.append(morph)
                elif "-interactions_distribution_metrics" in stem:
                    dist = pd.read_csv(file)
                    int_dist.append(dist)
                elif "-interactions_degree_metrics" in stem:
                    degree = pd.read_csv(file)
                    int_degree.append(degree)
                else:
                    print(f"File {stem} not recognized as interaction quantification data; skipping.")

    print(f"Found {fl_count} files from {ds_count} dataset(s) across {len(csv_path_list)} location(s).")

    # combine the org_morph lists found above into one combined table with all data
    labs_df = pd.concat(int_labs, axis=0, join='outer') if int_labs else None
    morph_df = pd.concat(int_morph, axis=0, join='outer') if int_morph else None
    dist_df = pd.concat(int_dist, axis=0, join='outer') if int_dist else None
    degree_df = pd.concat(int_degree, axis=0, join='outer') if int_degree else None

    # list all possible interaction site combinations
    all_pos = all_combos(organelle_names, splitter)

    ################################################
    # Summarize interactions count & morphology data
    ################################################
    if morph_df is not None:
        ### calculate interaction count/volume & summarize per organelle object for all interaction sites
        per_org_summary = perorg_interactions_cnt(interaction_morpho_df=morph_df, 
                                                   org_list=organelle_names,
                                                   splitter=splitter)

        # summarization parameters
        count_vol_group_by = ["dataset", "image_name", "mask_name", "scale", "object"]
        count_vol_cols = [col for col in per_org_summary.columns if col.endswith(("_count", "_volume"))]
        count_vol_ag_func_standard = {"num_interaction_types": ['mean', 'median', 'std']} | {col: ['sum', 'mean', 'median', 'std'] for col in count_vol_cols}

        # summarize per organelle type per image
        org_sum_tab = per_org_summary.groupby(count_vol_group_by).agg(count_vol_ag_func_standard)
    
        # Ensure all possible interactions are represented (if missing fill with NaN)
        for ind in org_sum_tab.index.droplevel(4).unique().to_list():
            for row in organelle_names:
                if ind+(row,) not in org_sum_tab.index:
                    org_sum_tab.loc[ind+(row,)] = np.nan
        org_sum_tab.sort_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-per_inter_count_volume_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-per_inter_count_volume_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            org_sum_tab.to_csv(str(out_path) + f"/{out_prefix}-per_inter_count_volume_summarystats.csv", mode='x')

        # unstack and format interaction count/volume summary table
        inter_count_vol_final = org_sum_tab.unstack(-1)
        for col in inter_count_vol_final.columns:
            if col[0].endswith(('_count', '_volume')):
                if col[2] not in col[0]:
                    inter_count_vol_final.drop(col,axis=1, inplace=True)

        inter_count_vol_final.columns = ["_".join((col_name[1], col_name[0], "per", col_name[-1])) for col_name in inter_count_vol_final.columns.to_flat_index()]
        inter_count_vol_final.columns = [col.replace('sum', 'total') for col in inter_count_vol_final.columns]
        inter_count_vol_final.columns = [col.replace('per', 'in') if 'total' in col else col for col in inter_count_vol_final.columns]
        inter_count_vol_final.fillna(0, inplace=True)
        inter_count_vol_final.reset_index()


        ### summarize interaction morphology per interaction site
        # summarization paramters
        morph_group_by = ["dataset", "image_name", "mask_name", "scale", "object"]
        morph_cols = ["SA_to_volume_ratio", "equivalent_diameter", "extent", "euler_number", "solidity", "axis_major_length"] + list(morph_df.filter(regex=".*intensity.*").columns)
        morph_ag_func_standard = ['mean', 'median', 'std']

        # summarize counts of interaction sites per image
        tab1 = morph_df[morph_group_by + ['ID']].groupby(morph_group_by).agg(['count'])
        tab1.rename(columns={'ID': 'sites'}, inplace=True)
        tab2 = morph_df.copy()[morph_df['in_higher_order'] == True][morph_group_by + ['ID']].groupby(morph_group_by).agg(['count'])
        tab2.rename(columns={'ID': 'sites_in_higher_order'}, inplace=True)
        tab3 = morph_df.copy()[morph_df['in_higher_order'] == False][morph_group_by + ['ID']].groupby(morph_group_by).agg(['count'])
        tab3.rename(columns={'ID': 'sites_not_in_higher_order'}, inplace=True)
        inter_sum_tab = pd.merge(tab1, tab2, 'outer', on=morph_group_by)
        inter_sum_tab = pd.merge(inter_sum_tab, tab3, 'outer', on=morph_group_by)

        # summarize all interaction sites
        tab4 = morph_df[morph_group_by + ['volume', 'surface_area']].groupby(morph_group_by).agg(['sum'] + morph_ag_func_standard)
        tab5 = morph_df[morph_group_by+morph_cols].groupby(morph_group_by).agg(morph_ag_func_standard)
        inter_sum_tab = pd.merge(inter_sum_tab, tab4, 'outer', on=morph_group_by)
        inter_sum_tab = pd.merge(inter_sum_tab, tab5, 'outer', on=morph_group_by)

        # Get mask_name and corresponding volume column per group & calculate volume fraction
        mask_names = morph_df.groupby(morph_group_by)['mask_name'].first()
        mask_volume_data = morph_df.groupby(morph_group_by).first().apply(lambda row: row[f"{mask_names.loc[row.name]}_volume"], axis=1)
        inter_sum_tab.insert(inter_sum_tab.columns.get_loc(('volume', 'sum')) + 1, ('volume', 'fraction'), inter_sum_tab[('volume', 'sum')]/mask_volume_data)

        # Ensure all possible interactions are represented (if missing fill with NaN)
        for ind in inter_sum_tab.index.droplevel(4).unique().to_list():
            for row in all_pos:
                if ind+(row,) not in inter_sum_tab.index:
                    inter_sum_tab.loc[ind+(row,)] = np.nan

        # fill NA with 0 for specific columns
        fill_dict = {('sites', 'count'): 0, 
                    ('sites_in_higher_order', 'count'): 0, 
                    ('sites_not_in_higher_order', 'count'): 0,
                    ('volume', 'sum'): 0,
                    ('surface_area', 'sum'): 0,
                    ('volume', 'fraction'): 0}
        inter_sum_tab = inter_sum_tab.fillna(value=fill_dict)

        # if (sites, count) is 1, set mean, median, and std to NaN
        single_site_mask = inter_sum_tab[('sites', 'count')] == 1
        for col in morph_cols+['volume', 'surface_area']:
            inter_sum_tab.loc[single_site_mask, (col, 'std')] = np.nan

        inter_sum_tab.sort_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-per_inter_morphology_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-per_inter_morphology_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            inter_sum_tab.to_csv(str(out_path) + f"/{out_prefix}-per_inter_morphology_summarystats.csv", mode='x')

        # unstack and format interaction morphology summary table   
        inter_morph_final = inter_sum_tab.unstack(-1)
        inter_morph_final.columns = ["_".join((col_name[1], col_name[-1], col_name[0])) for col_name in inter_morph_final.columns.to_flat_index()]
        inter_morph_final.columns = [col.replace('sum', 'total') for col in inter_morph_final.columns]
        inter_morph_final.columns = [col.replace(col, 'mask_volume') if 'mask' in col else col for col in inter_morph_final.columns]
        inter_morph_final = inter_morph_final.loc[:, ~inter_morph_final.columns.duplicated()]
        inter_morph_final.reset_index()

        # combine count/volume and morphology summaries
        final_combo_tab = pd.merge(inter_morph_final, inter_count_vol_final, on=["dataset", "image_name", "mask_name", "scale"]).reset_index()
    else:
        final_combo_tab = pd.DataFrame()

    ###################################
    # Summarize interaction labels data
    ###################################
    if labs_df is not None:
        ### summarize interaction site counts
        # define summarization paramters
        labs_group_by = ["dataset", "image_name", "mask_name", "object"]

        # summarize counts of interaction sites per image
        labs_tab1 = labs_df[labs_group_by + ['ID']].groupby(labs_group_by).agg(['count'])
        labs_tab1.rename(columns={'ID': 'sites'}, inplace=True)
        labs_tab2 = labs_df.copy()[labs_df['in_higher_order'] == True][labs_group_by + ['ID']].groupby(labs_group_by).agg(['count'])
        labs_tab2.rename(columns={'ID': 'sites_in_higher_order'}, inplace=True)
        labs_tab3 = labs_df.copy()[labs_df['in_higher_order'] == False][labs_group_by + ['ID']].groupby(labs_group_by).agg(['count'])
        labs_tab3.rename(columns={'ID': 'sites_not_in_higher_order'}, inplace=True)
        labs_inter_sum_tab = pd.merge(labs_tab1, labs_tab2, 'outer', on=labs_group_by)
        labs_inter_sum_tab = pd.merge(labs_inter_sum_tab, labs_tab3, 'outer', on=labs_group_by)

        # Ensure all possible interactions (all_pos) are represented (if missing fill with NaN):
        for ind in labs_inter_sum_tab.index.droplevel(3).unique().to_list():
            for row in all_pos:
                if ind+(row,) not in labs_inter_sum_tab.index:
                    labs_inter_sum_tab.loc[ind+(row,)] = np.nan

        # fill NA with 0 for specific columns
        fill_dict = {('sites', 'count'): 0, 
                    ('sites_in_higher_order', 'count'): 0, 
                    ('sites_not_in_higher_order', 'count'): 0}
        labs_inter_sum_tab = labs_inter_sum_tab.fillna(value=fill_dict)

        labs_inter_sum_tab.sort_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-per_inter_labels_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-per_inter_labels_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            labs_inter_sum_tab.to_csv(str(out_path) + f"/{out_prefix}-per_inter_labels_summarystats.csv", mode='x')

        # unstack and format interaction labels summary table
        inter_labels_final = labs_inter_sum_tab.unstack(-1)
        inter_labels_final.columns = ["_".join((col_name[1], col_name[-1], col_name[0])) for col_name in inter_labels_final.columns.to_flat_index()]
        inter_labels_final.reset_index()

        # combine with previous summary table
        final_combo_tab = pd.concat([final_combo_tab, inter_labels_final.reset_index()], axis=0)
    else:
        final_combo_tab = final_combo_tab                                                                            


    #########################################
    # Summarize interaction distribution data
    ########################################
    if dist_df is not None:
        ### summarize interaction site distribution metrics
        mask_name = "whole_image" if mask_name is None else mask_name

        if 'XY_center_vox_cnt_perbin' in list(dist_df.columns): # if there is a centering object
            nuc_dist_df = dist_df[["dataset", "image_name", 'mask_name', 'scale',
                                "XY_bins", "XY_center_vox_cnt_perbin", f"XY_{mask_name}_vox_cnt_perbin", "XY_center_cv_perbin",
                                "XY_wedges", "XY_center_vox_cnt_perwedge", f"XY_{mask_name}_vox_cnt_perwedge",
                                "Z_slices", "Z_center_vox_cnt", f"Z_{mask_name}_vox_cnt"]].drop_duplicates(subset=['dataset', 'image_name'])
            nuc_dist_df.columns = nuc_dist_df.columns.str.replace('center', 'obj', regex=False)
            nuc_dist_df.insert(loc=3,column='object',value='nuc')
            nuc_dist_df.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object'], inplace=True)


            inter_dist_df = dist_df[list(nuc_dist_df.reset_index().columns)]
            inter_dist_df.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object'], inplace=True)
            combo_dist_df = pd.concat([nuc_dist_df, inter_dist_df], axis=0)
        else: # if there is not a centering object
            dist_df.set_index(['dataset', 'image_name', "mask_name", 'scale', 'object'], inplace=True)
            combo_dist_df = dist_df

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

            dfs = [selection[['dataset', 'image_name', 'mask_name', 'scale', 'object']].reset_index()]
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
            for row in all_pos:
                if ind+(row,) not in dist_summary.index:
                    dist_summary.loc[ind+(row,)] = np.nan

        dist_summary.reset_index(inplace=True)

        # export before unstacking
        if (Path(out_path) / f"{out_prefix}-per_inter_distribution_summarystats.csv").exists():
            raise FileExistsError(f"CAUTION: {out_prefix}-per_inter_distribution_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
        else:
            dist_summary.to_csv(str(out_path) + f"/{out_prefix}-per_inter_distribution_summarystats.csv", mode='x')

        # unstack and format interaction distribution summary table
        dist_final = dist_summary.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object']).unstack(-1)
        dist_final.columns = ["_".join((col_name[1], col_name[0])) for col_name in dist_final.columns.to_flat_index()]
        dist_final = dist_final.reset_index()

        # combine with previous summary table
        final_combo_tab = pd.merge(final_combo_tab, dist_final, on=["dataset", "image_name", "mask_name", "scale"], how="outer")
    else:
        final_combo_tab = final_combo_tab

    ###################################
    # Summarize interaction degree data
    ###################################
    if degree_df is not None:
        # degree_df not exported before unstacking becuase it is already summarized per image originally

        # unstack and format interaction degree summary table
        inter_degree_final = degree_df.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object']).unstack(-1)
        inter_degree_final.columns = ["_".join((col_name[1], col_name[0])) for col_name in inter_degree_final.columns.to_flat_index()]
        inter_degree_final = inter_degree_final.reset_index()

        # combine with previous summary table
        final_combo_tab = pd.merge(final_combo_tab, inter_degree_final, on=["dataset", "image_name", "mask_name", "scale"], how="outer")
    else:
        final_combo_tab = final_combo_tab

    ##########################
    # Export combined results
    ##########################
    if (Path(out_path) / f"{out_prefix}-interactions_combined_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}-interactions_combined_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
    else:
        final_combo_tab.to_csv(str(out_path) + f"/{out_prefix}-interactions_combined_summarystats.csv", mode='x')

    print(f"Interactions summary is complete.")
    return final_combo_tab