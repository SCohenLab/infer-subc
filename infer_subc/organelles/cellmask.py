from typing import Dict
from pathlib import Path
import time
import numpy as np

from skimage.filters import scharr
from skimage.morphology import (isotropic_opening, 
                                isotropic_dilation, 
                                isotropic_erosion, 
                                binary_opening, 
                                binary_dilation, 
                                binary_erosion)

from skimage.morphology.footprints import ball, disk
from skimage.measure import label

from scipy.ndimage import zoom


from infer_subc.core.img import label_bool_as_uint16, make_aggregate, size_filter_linear_size
from infer_subc.core.file_io import export_inferred_organelle, import_inferred_organelle
from infer_subc.core.img import (
    masked_object_thresh,
    log_transform,
    min_max_intensity_normalization,
    scale_and_smooth,
    weighted_aggregate,
    masked_inverted_watershed,
    fill_and_filter_linear_size,
    get_max_label,
    get_interior_labels,
    select_cellmask_from_img
)

def raw_cellmask_fromaggr(img_in: np.ndarray, scale_min_max: bool = True) -> np.ndarray:
    """define cellmask image
    CELLMASK_W = (6.,1.,2.)
    CELLMASK_CH = (LYSO_CH,ER_CH,GOLGI_CH)

    Parameters
    ------------
    img_in
        a 3d image
    scale_min_max:
        scale to [0,1] if True. default True

    Returns
    -------------
        np.ndarray scaled aggregate

    """
    weights = (0, 6, 0, 2, 0, 1)
    if scale_min_max:
        return min_max_intensity_normalization(weighted_aggregate(img_in, *weights))
    else:
        return weighted_aggregate(img_in, *weights)

### USED ###
def non_linear_cellmask_transform(in_img):
    """non-linear distortion to fill out cellmask
    log + edge of smoothed composite

    Parameters
    ------------
    in_img:
        a 3d image containing all the channels

    Returns
    -------------
        np.ndarray scaled aggregate
    """
    # non-Linear processing
    log_img, d = log_transform(in_img.copy())
    log_img = min_max_intensity_normalization(log_img)
    return min_max_intensity_normalization(scharr(log_img)) + log_img

### USED ###
def choose_max_label_cellmask_union_nucleus(cellmask_img: np.ndarray, 
                                            cellmask_obj: np.ndarray, 
                                            nuclei_labels: np.ndarray, 
                                            watershed_method: str = 'slice-by-slice',
                                            interior_labels_only: bool = True
                                            ) -> np.ndarray:
    """get cellmask UNION nuclei for largest signal label

        Parameters
    ------------
    cellmask_img:
        the cellmask image intensities
    cellmask_obj:
        thresholded cellmask mask
    nuclei_labels:
        inferred nuclei labels (np.uint16)
    watershed_method:
        determines if the watershed should be run 'sice-by-slice' or in '3D' 

    Returns
    -------------
        boolean np.ndarray of cellmask+nuc corresponding to the label of largest total cellmask signal

    """

    cellmask_labels = masked_inverted_watershed(cellmask_img, nuclei_labels, cellmask_obj, method=watershed_method)

    # should we restrict to interior nuclear labels?
    # get_interior_labels(nuclei_object)
    # would need to update get_max_label to only choose the labels in get_interior_label
    target_labels = get_interior_labels(nuclei_labels) if interior_labels_only else None

    keep_label = get_max_label(cellmask_img, cellmask_labels, target_labels=target_labels)

    cellmask_out = np.zeros_like(cellmask_labels)
    cellmask_out[cellmask_labels == keep_label] = 1
    cellmask_out[nuclei_labels == keep_label] = 1

    return cellmask_out > 0

### USED ###
##########################
# 1. infer_cellmask
##########################
##########################
# infer_cellmask_fromaggr
##########################
def infer_cellmask_fromcomposite(in_img: np.ndarray,
                                  weights: list[int],
                                  rescale: bool,
                                  nuclei_labels: np.ndarray,
                                  median_sz: int,
                                  gauss_sig: float,
                                  mo_method: str,
                                  mo_adjust: float,
                                  mo_cutoff_size: int,
                                  min_hole_w: int,
                                  max_hole_w: int,
                                  small_obj_w: int,
                                  fill_filter_method: str,
                                  watershed_method: str
                                  ) -> np.ndarray:
    """
    Procedure to infer cellmask from linear unmixed input.

    Parameters
    ------------
    in_img: 
        a 3d image containing all the channels
    weights:
        a list of int that corresond to the weights for each channel in the composite; use 0 if a channel should not be included in the composite image
    rescale:
        True - rescale composite image
        False - don't rescale composite image
    nuclei_labels: 
        a 3d image containing the inferred nuclei labels
    median_sz: 
        width of median filter for _cellmask_ signal
    gauss_sig: 
        sigma for gaussian smoothing of _cellmask_ signal
    mo_method: 
         which method to use for calculating global threshold. Options include:
         "triangle" (or "tri"), "median" (or "med"), and "ave_tri_med" (or "ave").
         "ave" refers the average of "triangle" threshold and "mean" threshold.
    mo_adjust: 
        Masked Object threshold `local_adjust`
    mo_cutoff_size: 
        Masked Object threshold `size_min`
    max_hole_w: 
        hole filling cutoff for cellmask signal post-processing
    small_obj_w: 
        minimum object size cutoff for cellmask signal post-processing
    fill_filter_method:
        determines if the fill and filter function should be run 'sice-by-slice' or in '3D' 
    watershed_method:
        determines if the watershed should be run 'sice-by-slice' or in '3D' 

    Returns
    -------------
    cellmask_mask:
        a logical/labels object defining boundaries of cellmask

    """
    ###################
    # EXTRACT
    ###################
    struct_img = make_aggregate(in_img, *weights, rescale)

    ###################
    # PRE_PROCESSING
    ###################                         
    struct_img = scale_and_smooth(struct_img,
                                   median_size = median_sz, 
                                   gauss_sigma = gauss_sig)
    

    struct_img_non_lin = non_linear_cellmask_transform(struct_img)

    ###################
    # CORE_PROCESSING
    ###################
    struct_obj = masked_object_thresh(struct_img_non_lin, 
                                      global_method=mo_method, 
                                      cutoff_size=mo_cutoff_size, 
                                      local_adjust=mo_adjust)               

    ###################
    # POST_PROCESSING
    ###################
    struct_obj = fill_and_filter_linear_size(struct_obj, 
                                             hole_min=min_hole_w, 
                                             hole_max=max_hole_w, 
                                             min_size= small_obj_w,
                                             method=fill_filter_method)

    ###################
    # POST- POST_PROCESSING
    ###################
    cellmask_out = choose_max_label_cellmask_union_nucleus(struct_img, 
                                                           struct_obj, 
                                                           nuclei_labels, 
                                                           watershed_method=watershed_method) 

    return label_bool_as_uint16(cellmask_out)


# def fixed_infer_cellmask_fromcomposite(in_img: np.ndarray, nuclei_labels: np.ndarray) -> np.ndarray:
#     """
#     Procedure to infer cellmask from linearly unmixed input, with a *fixed* set of parameters for each step in the procedure.  i.e. "hard coded"

#     Parameters
#     ------------
#     in_img: 
#         a 3d image containing all the channels
#     nuclei_labels: 
#         a 3d image containing the inferred nuclei

#     Returns
#     -------------
#     cellmask_mask:
#         a logical/labels object defining boundaries of cellmask
#     """
    

#     ###################
#     # PARAMETERS
#     ###################   
#     weights = [0,0,0,3,3,2]
#     median_sz = 10
#     gauss_sig = 1.34
#     mo_method = "med"
#     mo_adjust = 0.3
#     mo_cutoff_size = 150
#     hole_min_width = 0
#     hole_max_width = 50
#     small_obj_w = 45
#     fill_filter_method = '3D'
#     watershed_method = '3D'

#     cellmask_out = infer_cellmask_fromcomposite(in_img,
#                                                 weights,
#                                                 nuclei_labels,
#                                                 median_sz,
#                                                 gauss_sig,
#                                                 mo_method,
#                                                 mo_adjust,
#                                                 mo_cutoff_size,
#                                                 hole_min_width,
#                                                 hole_max_width,
#                                                 small_obj_w,
#                                                 fill_filter_method,
#                                                 watershed_method) 
    
#     return cellmask_out.astype(np.uint8)


##########################
#  infer_cellmask_fromcytoplasm
##########################
def infer_cellmask_fromcytoplasm(cytoplasm_mask: np.ndarray,
                                  nucleus_mask: np.ndarray,
                                  min_hole_width: int,
                                  max_hole_width: int,
                                  small_obj_width: int,
                                  fill_filter_method: str
                                  ) -> np.ndarray:
    """
    Procedure to infer 3D nuclei segmentation from multichannel z-stack input.

    Parameters
    ------------
    cytoplasm_mask: np.ndarray
        3D image containing the mask of the cytoplasm
    nucleus_mask: np.ndarray
        3D image containing the mask of the nucleus
    min_hole_width: int
        minimum size of holes to fill in final cell mask
    max_hole_width: int,
        maximum size of holes to fill in final cell mask
    small_obj_w: int
        minimum object size cutoff to remove from final cell mask; likely not required since small objects were removed from cytoplasm mask
    fill_method: str
        method for fill and filter; either "3D" or "slice_by_slice"

    Returns
    -------------
    cell_mask
        mask defined extent of the entire cell
    
    """

    ###################
    # CORE_PROCESSING
    ###################
    cell = np.logical_or(nucleus_mask, cytoplasm_mask)

    ###################
    # POST_PROCESSING
    ###################
    cleaned_img = fill_and_filter_linear_size(cell, 
                                              hole_min=min_hole_width, 
                                              hole_max=max_hole_width, 
                                              min_size=small_obj_width, 
                                              method=fill_filter_method)

    ###################
    # LABELING
    ###################
    cell_mask = label_bool_as_uint16(cleaned_img)

    return cell_mask


##########################
#  fixed_infer_cellmask_fromcytoplasm
##########################
# def fixed_infer_cellmask_fromcytoplasm(cytoplasm_mask: np.ndarray,
#                                         nucleus_mask:np.ndarray) -> np.ndarray:
#     """
#     Procedure to infer cellmask from the cytoplasm mask

#     Parameters
#     ------------
#     in_img: np.ndarray
#         a 3d image containing cytoplasm segmentation
 
#     Returns
#     -------------
#     nuclei_object
#         inferred nuclei
    
#     """
#     min_hole_w = 0
#     max_hole_w = 30
#     small_obj_w = 0
#     fill_filter_method = "3D"

#     return infer_cellmask_fromcytoplasm(cytoplasm_mask,
#                                          nucleus_mask,
#                                          min_hole_w,
#                                          max_hole_w,
#                                          small_obj_w,
#                                          fill_filter_method)




# def infer_and_export_cellmask(
#     in_img: np.ndarray, nuclei_obj: np.ndarray, meta_dict: Dict, out_data_path: Path
# ) -> np.ndarray:
#     """
#     infer cellmask and write inferred cellmask to ome.tif file

#     Parameters
#     ------------
#     in_img:
#         a 3d  np.ndarray image of the inferred organelle (labels or boolean)
#     nuclei_obj:
#         a 3d image containing the inferred nuclei
#     meta_dict:
#         dictionary of meta-data (ome)
#     out_data_path:
#         Path object where tiffs are written to

#     Returns
#     -------------
#     exported file name

#     """
#     cellmask = fixed_infer_cellmask_fromcomposite(in_img, nuclei_obj)
#     out_file_n = export_inferred_organelle(cellmask, "cell", meta_dict, out_data_path)
#     print(f"inferred cellmask. wrote {out_file_n}")
#     return cellmask>0


# def get_cellmask(in_img: np.ndarray, nuclei_obj: np.ndarray, meta_dict: Dict, out_data_path: Path) -> np.ndarray:
#     """
#     load cellmask if it exists, otherwise calculate and write inferred cellmask to ome.tif file

#     Parameters
#     ------------
#     in_img:
#         a 3d  np.ndarray image of the inferred organelle (labels or boolean)
#     nuclei_obj:
#         a 3d image containing the inferred nuclei
#     meta_dict:
#         dictionary of meta-data (ome)
#     out_data_path:
#         Path object where tiffs are written to

#     Returns
#     -------------
#     exported file name

#     """

#     try:
#         cellmask = import_inferred_organelle("cell", meta_dict, out_data_path)
#     except:
#         start = time.time()
#         print("starting segmentation...")
#         cellmask = fixed_infer_cellmask_fromcomposite(in_img, nuclei_obj)
#         out_file_n = export_inferred_organelle(cellmask, "cell", meta_dict, out_data_path)
#         end = time.time()
#         print(f"inferred (and exported) cellmask in ({(end - start):0.2f}) sec")

    # return cellmask

### USED ###
##########################
#  create the cell mask by adding the nuclei and cytoplasm masks together
##########################
def combine_cytoplasm_and_nuclei(cyto_seg: np.ndarray,
                                 nuc_seg: np.ndarray,
                                 max_hole_width: int):
    """
    Function to combine the the cytoplasm and nuclei segmentations to produce the entire cell mask.

    Parameters:
    ----------
    cyto_seg: np.ndarray,
        image containing the cytoplasm segmentation
    nuc_seg: np.ndarray,
        image containing the nuclei segmentation
    max_hole_width: int
        size of the gaps between the nuclei and cytoplasm (usually small)
    """ 
    
    cells = np.logical_or(cyto_seg.astype(bool), nuc_seg.astype(bool))

    cell_multiple = fill_and_filter_linear_size(cells, 
                                                hole_min=0,
                                                hole_max=max_hole_width,
                                                min_size=0,
                                                method='3D')
    
    cell_area = cell_multiple.astype(bool)

    return cell_area

### USED ###
def select_highest_intensity_cell(raw_image: np.ndarray,
                                   cell_seg: np.ndarray,
                                   nuc_seg: np.ndarray):
    """ 
    Create an instance segmentation of the cell area using a watershed operation based on nuclei seeds.
    Then, select the cell with the highest combined organelle intensity.

    Parameters:
    ----------
    raw_image: np.ndarray,
        gray scale 3D multi-channel numpy array (CZYX)
    cell_seg: np.ndarray,
        binary cell segmentation with multiple cells in the FOV
    nuc_seg: np.ndarray,
        labeled nuclei segmentation with each nuclei having a different ID number (e.g., the result of the skimage label() function)
    labels_to_consider: Union(list, None)
        a list of labels that should be considered when determining the highest intensity. Default is None which utilizes all possible labels in the cell image
        
    Output
    ----------
    good_cell: np.ndarray  
        a binary image of the single cell with the highest total fluorescence intensity
    """
    # instance segmentation of cell area with watershed function
    cell_labels = masked_inverted_watershed(cell_seg, markers=nuc_seg, mask=cell_seg, method='3D')

    # create composite of all fluorescence channels after min-max normalization
    norm_channels = [(min_max_intensity_normalization(raw_image[c])) for c in range(len(raw_image))]
    normed_signal = np.stack(norm_channels, axis=0)
    normed_composite = normed_signal.sum(axis=0)

    # list of cell IDs to measure intensity of
    all_labels = np.unique(cell_labels)[1:]

    # measure total intensity in each cell from the ID list
    total_signal = [normed_composite[cell_labels == label].sum() for label in all_labels]

    # select the cell with the highest total intensity
    keep_label = all_labels[np.argmax(total_signal)]
    good_cell = cell_labels == keep_label

    return good_cell

def find_radius(cell_mask: np.ndarray, method: str, verbose: bool=False) -> np.ndarray:
    """
    Determines the radius of the cells in the mask, and outputs a copy of the mask with the radii encoded as their labels.

    Parameters:
    ----------
    cell_mask : np.ndarray
        A mask of the cells.
    method : str
        The method to use for finding the radius. Can be 'isotropic' or 'binary'.
        The `binary` method uses a saucer-like shape to remove the soma from the neurites. 
        The "saucer" is comprised of 3 2-Dimensional disks stacked into a 3-Dimensional space. 
        This method is more likely to oversegment the soma as compared to the `isotropic` methods. 
        This results in the initial portions of the neurites possibly being included in the soma segmentation. 
        Using the `binary` method will allow for the selection of more irregularly shaped soma objects.
    verbose: bool
        If True, prints out information about the radius finding process.
        
    Returns:
    -------
    np.ndarray
        A mask of the cells with their radii encoded as labels.
    """
    radii_mask = np.zeros_like(cell_mask)
    cell_mask_resize = zoom(cell_mask.copy(), (1, 0.5, 0.5)) # resizing cell mask to speed up processing
    zz, yy, xx = cell_mask_resize.shape                      # collecting y-length of image NOTE: may want to collect sqrt(yy^2 + xx^2) and use in place of yy   

    cell_nums = np.unique(cell_mask[cell_mask != 0])
    label_factor = 10 ** len(str(cell_nums.max()))

    for cell_num in cell_nums:
        test_img = (cell_mask_resize == cell_num)
        rad_range = [i+1 for i in range(yy // 4)]  # Dividing by 4 because mask is resized

        if method == 'isotropic':
            while len(rad_range) > 2:                   # repeats code until only 1 or 2 radii remain
                rad = rad_range[len(rad_range) // 2]    # sets test radius to radius in middle of rad_range list
                if verbose: print(f"Trying radius of {rad}")

                # testing erosion with test radius
                if np.all(isotropic_erosion(test_img.astype(np.uint8), rad) == 0):
                    rad_range = rad_range[:rad_range.index(rad)]
                    if verbose: print(f"{rad} is too large")
                else:
                    rad_range = rad_range[rad_range.index(rad)+1:]
                    if verbose: print(f"{rad} is too small")
                if verbose: print(f"{len(rad_range)} possible radii remaining")
        elif method == 'binary':
            while len(rad_range) > 2:                   # repeats code until only 1 or 2 radii remain
                rad = rad_range[len(rad_range) // 2]    # sets test radius to radius in middle of rad_range list
                if verbose: print(f"Trying radius of {rad}")

                #creates 'saucer'
                edge = disk(rad // 4)
                middle = disk(rad)
                w = (middle.shape[0] - edge.shape[0]) // 2
                edge = np.pad(edge, ((w, w), (w, w)), mode='constant', constant_values=0)
                fp = np.stack((edge, middle, edge))

                # testing erosion using test radius
                if np.all(binary_erosion(test_img.astype(np.uint8), fp) == 0):
                    rad_range = rad_range[:rad_range.index(rad)]
                    if verbose: print(f"{rad} is too large")
                else:
                    rad_range = rad_range[rad_range.index(rad)+1:]
                    if verbose: print(f"{rad} is too small")
                if verbose: print(f"{len(rad_range)} possible radii remaining")

        if len(rad_range) == 1:
            opti_rad = rad_range[0] // 2
        elif len(rad_range) == 2:
            opti_rad = (rad_range[0] + rad_range[1]) // 4
        radii_mask[cell_mask == cell_num] = (opti_rad * label_factor) + cell_num
    return radii_mask

def infer_soma_from_mask(cell_mask: np.ndarray, radii_mask: np.ndarray, method: str='binary'):
    """
    Infers the soma region from the cell mask and radii mask by deriving the radius of each cell.

    Parameters:
    ----------

    cell_mask : np.ndarray
        A mask of the cells.
    radii_mask : np.ndarray
        A mask of the cells with their radii encoded as labels.
    method : str
        The method to use for inferring the soma. Can be 'isotropic' or 'binary'.

    Returns:
    -------
    np.ndarray
        A mask of the inferred soma regions.
    """
    soma_out_1 = np.zeros_like(cell_mask)

    cell_nums = np.unique(cell_mask[cell_mask != 0])
    label_factor = 10 ** len(str(cell_nums.max()))

    for cell_num in cell_nums:                              # repeat for each cell in the image
        # determine the radius of the chosen cell
        soma_img_solo = (cell_mask == cell_num)
        opti_rad = np.unique(radii_mask[soma_img_solo])[0]
        opti_rad = (opti_rad - cell_num) / label_factor

        if method == 'isotropic': 
            # opening and dilation to ensure removal of neurites
            neurites_removed = isotropic_opening(soma_img_solo.astype(np.uint8), opti_rad)
            soma_initial = isotropic_dilation(neurites_removed, opti_rad) & soma_img_solo
        elif method == 'binary':

            # creation of 'saucer'
            edge = disk(int(opti_rad // 2))
            middle = disk(int(opti_rad))
            w = (middle.shape[0] - edge.shape[0]) // 2
            edge = np.pad(edge, ((w, w), (w, w)), mode='constant', constant_values=0)
            fp = np.stack((edge, middle, edge))

            #opening and dilation to ensure removal of neurites
            neurites_removed = binary_opening(soma_img_solo.astype(np.uint8), fp)
            soma_initial = binary_dilation(neurites_removed, footprint=ball(int(opti_rad // 2))) & soma_img_solo
        else:
            raise ValueError(f"method of {method} was given, but only 'isotropic' or 'binary' is allowed.")

        soma_out_1[soma_initial] = cell_num
    return soma_out_1

def infer_neurites_from_mask(cell_mask: np.ndarray, radii_mask: np.ndarray, soma_out_1: np.ndarray, method: str):
    """
    Uses the cell mask and the soma output to infer the neurite regions.

    Parameters:
    ----------

    cell_mask : np.ndarray
        A mask of the cells.
    radii_mask : np.ndarray
        A mask of the cells with their radii encoded as labels.
    soma_out_1 : np.ndarray
        A mask of the inferred soma regions.
    method : str
        The method to use for inferring the neurites. Can be 'isotropic' or 'binary'.

    Returns:
    -------
    np.ndarray
        A mask of the inferred neurite regions.
    """
    neurites_out_1 = np.zeros_like(cell_mask)

    cell_nums = np.unique(cell_mask[cell_mask != 0])
    label_factor = 10 ** len(str(cell_nums.max()))
    binary_soma = soma_out_1 > 0

    for cell_num in cell_nums:                          # repeat across cell numbers
        # determine radius for cell
        solo_mask = (cell_mask == cell_num)
        opti_rad = np.unique(radii_mask[solo_mask])[0]
        opti_rad = (opti_rad - cell_num) / label_factor

        # mask out soma from full cell mask
        neurite_mask = ~binary_soma & solo_mask

        # filter out small objects that may instead be missing outcrops from the soma
        # size exclusion was determined imperically on example neurons from 63X magnification images.
        # TODO: add optional parameter to adjust size filteringre
        if method == 'isotropic':
            filtered = size_filter_linear_size(img=label(neurite_mask), min_size=(opti_rad*2), method='3D') * solo_mask
        elif method == 'binary':
            filtered = size_filter_linear_size(img=label(neurite_mask), min_size=(opti_rad//2), method='3D') * solo_mask
        else:
            raise ValueError(f"method of {method} was given, but only 'isotropic' or 'binary' is allowed.")

        # label the neurites to unique IDs while also encoding their cell radius
        neurite_labels = label(filtered)
        neurite_labels[neurite_labels > 0] = (neurite_labels[neurite_labels > 0] * label_factor) + cell_num
        neurites_out_1[solo_mask] = neurite_labels[solo_mask]
    return neurites_out_1

def clean_soma_from_neurites(cell_mask: np.ndarray, neurites_out_1: np.ndarray) -> np.ndarray:
    """
    Cleans the soma regions from the neurites by masking out the neurites from the cell mask.

    Parameters:
    ----------

    cell_mask : np.ndarray
        A mask of the cells.
    neurites_out_1 : np.ndarray
        A mask of the inferred neurite regions.

    Returns:
    -------
    np.ndarray
        A cleaned mask of the soma regions.
    """
    soma_out_2 = np.zeros_like(cell_mask)

    cell_nums = np.unique(cell_mask[cell_mask != 0])
    label_factor = 10 ** len(str(cell_nums.max()))

    # Create a mask for all neurites at once
    neurites_mask = (neurites_out_1 % label_factor) > 0

    # For each cell, mask soma regions in one go
    soma_mask = (~neurites_mask) & (cell_mask != 0)

    # Find the most common value in soma_mask for each cell and assign only those pixels
    for cell_num in cell_nums:
        cell_region = (cell_mask == cell_num)
        soma_region = label(soma_mask & cell_region)
        # Only keep the largest connected region (most common value)
        if np.any(soma_region):
            bincount = np.bincount(soma_region.ravel())
            main_val = np.argmax(bincount[1:]) + 1 if len(bincount) > 1 else 1
            soma_region = soma_region & (soma_region == main_val)
            soma_out_2[cell_region] = soma_region[cell_region] * cell_num
    return soma_out_2

def clean_neurites_from_soma(cell_mask: np.ndarray, soma_out_2: np.ndarray):
    """
    Cleans the neurite regions from the soma by masking out the soma from the neurite mask.

    Parameters:
    ----------
    cell_mask : np.ndarray
        A mask of the cells.
    soma_out_2 : np.ndarray
        A mask of the cleaned soma regions.

    Returns:
    -------
    np.ndarray
        A cleaned mask of the neurite regions.
    """
    neurites_out_2 = np.zeros_like(cell_mask)

    cell_nums = np.unique(cell_mask[cell_mask != 0])
    binary_soma = soma_out_2 > 0

    # Create a mask for all neurites at once
    neurites_mask = (~binary_soma) & (cell_mask != 0)

    # Label all neurite regions in one call
    neurites_labels = label(neurites_mask)

    # Relabel to encode cell number
    label_factor = 10 ** len(str(cell_nums.max()))
    neurites_out_2 = np.zeros_like(cell_mask)
    for cell_num in cell_nums:
        neurites_out_2[(cell_mask == cell_num) & (neurites_mask > 0)] = (neurites_labels[(cell_mask == cell_num) & (neurites_mask > 0)] * label_factor) + cell_num
    return neurites_out_2

def infer_soma_neurites(in_seg: np.ndarray, multichannel_input: bool=False, chan: int=0, rad_method: str='method', soma_method: str='method', neurite_method: str='method', method='binary'):
    """
    Infers the soma and neurite regions from the input segmentation based on either the binary or isotropic filtering method.

    Parameters:
    ----------
    in_seg : np.ndarray
        The input segmentation mask.
    multichannel_input : bool
        Whether the input is a multichannel image.
    chan : int
        The channel to use for segmentation if it is a multichannel image.
    rad_method: str
        The method used to approximate the radius ('binary' or 'isotropic').
    soma_method: str
        The method used to approximate the soma ('binary' or 'isotropic').
    neurite_method: str
        The method used to approximate the neurites ('binary' or 'isotropic').
    method : str
        The optional method to used for all unassigned methods ('binary' or 'isotropic') (optional).

    Returns:
    -------
    np.ndarray
        A stack of the inferred soma and neurite regions.
    """
    ###################
    # EXTRACT
    ###################  

    if (rad_method == 'method'):
        rad_method = method
    if (soma_method == 'method'):
        soma_method = method
    if (neurite_method == 'method'):
        neurite_method = method

    cell_mask = select_cellmask_from_img(in_seg, multichannel_input=multichannel_input, chan=chan)

    ###################
    # PRE_PROCESSING
    ################### 
    radii_mask = find_radius(cell_mask, rad_method)

    ###################
    # CORE_PROCESSING
    ###################
    soma_initial = infer_soma_from_mask(cell_mask, radii_mask, soma_method)

    neurites_initial = infer_neurites_from_mask(cell_mask, radii_mask, soma_initial, neurite_method)

    ###################
    # POST_PROCESSING
    ################### 
    soma_cleaned = clean_soma_from_neurites(cell_mask, neurites_initial)

    neurites_cleaned = clean_neurites_from_soma(cell_mask, soma_cleaned)

    ###################
    # POST_POST_PROCESSING
    ################### 
    soma_neurites = np.stack([soma_cleaned, neurites_cleaned])
    
    return soma_neurites