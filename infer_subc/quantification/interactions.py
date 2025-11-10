
import pandas as pd
import numpy as np

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
    
    '''
    Determine which interaction sites are included in higher order interactions.
    An interaction site is considered to be part of a higher order interaction if it overlaps with an additional organelle not included in the original interaction site definition.
    For example, if the interaction site you are creating is mitoXlyso, and a specific interaction site in the image also overlaps with ER, then that interaction site is considered to be part of a higher order interaction.
    The output is an image array of the interaction sites with unique integer IDs for each site that is NOT part of a higher order interaction.

    Parameters
    ----------
    site : np.ndarray
        An image array of the overlap regions between the organelles included in the organelle_segs variable;
        each interaction site should be labeled with unique integer IDs that will be included in the output image.
    site_name : str
        A string of organelle names separated by the specified splitter.
    inter_labels_tab : pd.DataFrame
        A pandas DataFrame table with unique identifiers (integer IDs and labels) associated to each interaction site.
        `ID`: unique integer identifier for each interaction site in the overlap image. Each site will have a different ID number.
        `object`: the name of the interaction sites being examined, created by joining the organelle names with the specified splitter.
        `label`: a string of organelle ID numbers involved in each interaction site, joined by underscores. 
    organelle_segs : dict[str:np.ndarray]
        A dictionary of organelle segmentations for all organelles from the same cell, including the ones in the interaction site and other organelles to check against for higher order interactions.
        The dictionary has organelle names as keys and segmentation image arrays as values.
    splitter : str, optional
        The character used to split the organelle names in the orgs string, by default "X". 
        For example, "mitoXlyso" would indicate an interaction between mito and lyso.
    
    Returns
    -------
    LOi_NR : np.ndarray
        An image array of the interaction sites with unique integer IDs for each site that is NOT part of a higher order interaction.
    inter_labels_tab : pd.DataFrame
        The input pandas DataFrame table with an additional column indicating if the interaction site is part of a higher order interaction.
        `in_higher_order`: a boolean value where True/1 indicates the interaction site is part of a higher order interaction, and False/0 indicates it is not.
    '''
    # remove any interaction sites that are involved in higher order interactions
    LOc_NR = site.copy()            
    for org, val in organelle_segs.items():         
        if (org not in site_name.split(splitter)
            and np.any(site.astype(int)*val.astype(int))):
            HOc = site.copy()       
            valid = (LOc_NR>0)*(val>0)                  
            HOc[valid.astype(bool)==False]=0
            for id in np.unique(HOc):
                LOc_NR[LOc_NR==id] = 0    

    # ensure the original site IDs are preserved
    LOi_NR = (LOc_NR>0).astype(int) * site

    # select only the positive integer values within the array
    redundancy = inter_labels_tab['ID'].isin(np.unique(LOi_NR[LOi_NR>0]).tolist())

    # add new column to the interaction table indicating if the site is in a higher order interaction
    new_tab = inter_labels_tab.copy()
    new_tab.insert((inter_labels_tab.columns.get_loc('label')+1), "in_higher_order", list(map(bool, ~redundancy)))

    return LOi_NR, new_tab


def create_interaction_sites(org_name_list:List[str],
                             org_seg_list: List[np.ndarray],
                             name_splitter: str="X") -> tuple[np.ndarray, pd.DataFrame]:
    
    '''
    Create an image of the overlap regions between the selected organelles and a table of unique identifiers associated to each interaction site.
    
    Parameters 
    ----------
    org_name_list : List[str]
        A list of organelle names as strings. These will be the organelles used to create interaction sites.
    org_seg_list : List[np.ndarray]
        A list of organelle segmentation images as numpy ndarrays. These should be in the same order as the org_name_list list.
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
    
    if len(org_name_list)<2:
        raise ValueError("Please select at least two organelles to define an interaction site.")
    else:
        # create name for interaction site
        interaction_name = name_splitter.join(org_name_list)

        # run function to create dictionary of organelle segmentations
        org_dict = make_dict(org_name_list, org_seg_list)

        # create the overlap image
        overlap_img = create_overlap(interaction_name, org_dict)

        # use regionprops table to list interaction sites by unique index and extract slice for each object
        inter_tab = find_inter_labels(overlap_img, interaction_name, org_dict)

        # determine if each site is also involved in a higher order interaction (there are more than the specified organelles involved)
        lower_order_sites, inter_tab = assess_if_higher_order_int(overlap_img, interaction_name, inter_tab, org_dict)

        return overlap_img, lower_order_sites, inter_tab