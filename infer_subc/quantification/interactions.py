
import pandas as pd
import numpy as np
import math

from skimage.measure import regionprops_table

from infer_subc.core.img import *
from infer_subc.quantification.stats import *
from infer_subc.quantification.stats_helpers import *
from infer_subc.organelles import * 


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
                       org_dict: dict[str, np.ndarray]) -> pd.DataFrame:
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
    involved = interaction_name.split("X")
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
                              name_splitter: str="X") -> tuple[np.ndarray, pd.DataFrame]:
    
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
        A list of all organelle segmentation images as numpy ndarrays. These should be in the same order as the org_name_list list.
    name_splitter : str, optional
        The character used to separate the organelle names in the org_name_list string, by default "X". 
        For example, "mitoXlyso" would indicate an interaction between mito and lyso. 
        Use of other splitters may cause issues during downstream analysis in infer-subc. 
        Specifically, avoid using "_" or "-" as a splitter as they are used in other parts of the analysis.
    
    Returns
    -------
    overlap_img : np.ndarray
        An image array of the overlap regions between the selected organelles, with unique integer IDs for each interaction site.   
    inter_tab : pd.DataFrame
        A pandas DataFrame table with unique identifiers (integer IDs and labels) associated to each interaction site.
        `ID`: unique integer identifier for each interaction site in the overlap image. Each site will have a different ID number.
        `object`: the name of the interaction sites being examined, created by joining the organelle names with the specified splitter.
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

        # use regionprops table to list interaction sites by unique index and extract slice for each object
        inter_tab = find_inter_labels(overlap_img, interaction_name, org_dict)

        # determine if each site is also involved in a higher order interaction (there are more than the specified organelles involved)
        lower_order_sites, inter_tab = assess_if_higher_order_int(overlap_img, interaction_name, inter_tab, org_dict)

        return overlap_img, lower_order_sites, inter_tab
    

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