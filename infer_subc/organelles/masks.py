from typing import Union

import numpy as np

from infer_subc.core.img import *
from infer_subc.organelles.cellmask import (infer_cellmask_fromcomposite,
                                            infer_cellmask_fromcytoplasm,
                                            select_highest_intensity_cell,
                                            combine_cytoplasm_and_nuclei)

from infer_subc.organelles.nuclei import infer_nuclei_fromlabel, infer_nuclei_fromcytoplasm, mask_cytoplasm_nuclei, segment_nuclei_seeds
from infer_subc.organelles.cytoplasm import infer_cytoplasm, infer_cytoplasm_fromcomposite, infer_cytoplasm_fromcomposite
from infer_subc.organelles.membrane import infer_intermediate_masks, infer_cellmask_masks_C, infer_nucleus_masks_C, infer_cellmask_masks_D

##################
# Masks Workflow #
##################
def infer_masks(in_img: np.ndarray, 
                nuc_ch: Union[int,None],
                nuc_median_sz: int, 
                nuc_gauss_sig: float,
                nuc_thresh_factor: float,
                nuc_thresh_min: float,
                nuc_thresh_max: float,
                nuc_min_hole_w: int,
                nuc_max_hole_w: int,
                nuc_small_obj_w: int,
                nuc_fill_filter_method: str,
                cell_weights: list[int],
                cell_rescale: bool,
                cell_median_sz: int,
                cell_gauss_sig: float,
                cell_mo_method: str,
                cell_mo_adjust: float,
                cell_mo_cutoff_size: int,
                cell_min_hole_w: int,
                cell_max_hole_w: int,
                cell_small_obj_w: int,
                cell_fill_filter_method: str,
                cell_watershed_method: str,
                cyto_erode_nuclei = True
                
                ):
    """
    Procedure to infer nucleus, cellmask, and cytoplams from multichannel confocal microscopy input image.

    Parameters
    ------------
    in_img: np.ndarray
        a 3d image containing all the organelle and nuclei label channels
    nuc_median_sz: int
        width of median filter for nuclei channel
    nuc_gauss_sig: float
        sigma for gaussian smoothing for nuclei channel
    nuc_thresh_factor: float
        adjustment factor for log Li threholding for nuclei channel
    nuc_thresh_min: float
        abs min threhold for log Li threholding for nuclei channel
    nuc_thresh_max: float
        abs max threhold for log Li threholding for nuclei channel
    nuc_max_hole_w: int
        hole filling cutoff for nuclei post-processing for nuclei channel
    nuc_small_obj_w: int
        minimum object size cutoff for nuclei post-processing for nuclei channel
    cell_weights:
        a list of int that corresond to the weights for each channel in the composite; use 0 if a channel should not be included in the composite image
    cell_rescale:
        True - rescale composite image
        False - don't rescale composite image
    cell_nuclei_labels: 
        a 3d image containing the inferred nuclei labels
    cell_median_sz: 
        width of median filter for _cellmask_ signal
    cell_gauss_sig: 
        sigma for gaussian smoothing of _cellmask_ signal
    cell_mo_method: 
         which method to use for calculating global threshold. Options include:
         "triangle" (or "tri"), "median" (or "med"), and "ave_tri_med" (or "ave").
         "ave" refers the average of "triangle" threshold and "mean" threshold.
    cell_mo_adjust: 
        Masked Object threshold `local_adjust`
    cell_mo_cutoff_size: 
        Masked Object threshold `size_min`
    cell_max_hole_w: 
        hole filling cutoff for cellmask signal post-processing
    cell_small_obj_w: 
        minimum object size cutoff for cellmask signal post-processing
    cell_fill_filter_method:
        determines if the fill and filter function should be run 'sice-by-slice' or in '3D' 
    cell_watershed_method:
        determines if the watershed should be run 'sice-by-slice' or in '3D' 
    cyto_erode_nuclei: 
        should we erode? Default False
    
    Output
    ------
    stacked_masks: np.ndarray
        a multi-channel np.ndarry containing the nuclei, cell, and cytoplasm masks
    """

    ####################
    ### infer nuclei ###
    ####################
    nuc_objs = infer_nuclei_fromlabel(in_img,
                                      nuc_ch,
                                      nuc_median_sz,
                                      nuc_gauss_sig,
                                      nuc_thresh_factor,
                                      nuc_thresh_min, 
                                      nuc_thresh_max, 
                                      nuc_min_hole_w, 
                                      nuc_max_hole_w,
                                      nuc_small_obj_w, 
                                      nuc_fill_filter_method)
    
    ######################
    ### infer cellmask ###
    ######################
    cellmask_obj = infer_cellmask_fromcomposite(in_img,
                                                cell_weights,
                                                cell_rescale,
                                                nuc_objs,
                                                cell_median_sz,
                                                cell_gauss_sig,
                                                cell_mo_method,
                                                cell_mo_adjust,
                                                cell_mo_cutoff_size,
                                                cell_min_hole_w,
                                                cell_max_hole_w,
                                                cell_small_obj_w,
                                                cell_fill_filter_method,
                                                cell_watershed_method)

    ######################
    ### infer cellmask ###
    ###################### 
    cyto_obj = infer_cytoplasm(nuc_objs, 
                               cellmask_obj,
                               erode_nuclei=cyto_erode_nuclei)
    
    #########################
    ### select single nuc ###
    #########################
    nuc_obj = apply_mask(nuc_objs, cellmask_obj)
    nuc_obj = label_bool_as_uint16(nuc_obj)
    
    ###################
    ### stack masks ###
    ###################
    maskstack = stack_masks(nuc_mask=nuc_obj, 
                            cellmask=cellmask_obj,
                            cyto_mask=cyto_obj)
    
    return maskstack

####################
# Masks A Workflow #
####################
def infer_masks_A(in_img: np.ndarray,
                    cyto_weights: list[int],
                    cyto_rescale: bool,
                    cyto_median_sz: int,
                    cyto_gauss_sig: float,
                    cyto_mo_method: str,
                    cyto_mo_adjust: float,
                    cyto_mo_cutoff_size: int,
                    cyto_min_hole_w: int,
                    cyto_max_hole_w: int,
                    cyto_small_obj_w: int,
                    cyto_fill_filter_method: str,
                    nuc_min_hole_w: int,
                    nuc_max_hole_w: int,
                    nuc_fill_method: str,
                    nuc_small_obj_w: int,
                    nuc_fill_filter_method: str,
                    cell_min_hole_width: int,
                    cell_max_hole_width: int,
                    cell_small_obj_width: int,
                    cell_fill_filter_method: str
                    ) -> np.ndarray:
    """
    Procedure to infer cellmask from linear unmixed input.

    Parameters
    ------------
    in_img: 
        a 3d image containing all the channels
    cyto_weights:
        a list of int that corresond to the weights for each channel in the composite; use 0 if a channel should not be included in the composite image
    cyto_rescale:
        True = rescale composite
        False = don't rescale composite
    cyto_median_sz: 
        width of median filter for _cellmask_ signal
    cyto_gauss_sig: 
        sigma for gaussian smoothing of _cellmask_ signal
    cyto_mo_method: 
         which method to use for calculating global threshold. Options include:
         "triangle" (or "tri"), "median" (or "med"), and "ave_tri_med" (or "ave").
         "ave" refers the average of "triangle" threshold and "mean" threshold.
    cyto_mo_adjust: 
        Masked Object threshold `local_adjust`
    cyto_mo_cutoff_size: 
        Masked Object threshold `size_min`
    cyto_max_hole_w: 
        hole filling cutoff for cellmask signal post-processing
    cyto_small_obj_w: 
        minimum object size cutoff for cellmask signal post-processing
    cyto_fill_filter_method:
        determines if fill and filter should be run 'sice-by-slice' or in '3D' 
    nuc_max_hole_w: int
        hole filling cutoff to fill the nuclei
    nuc_small_obj_w: int
        object size cutoff to remove artifacts from dilation/erosion steps
    nuc_fill_filter_method: str
        to filter artifacts in "3D" or "slice-by-slice"
    cell_min_hole_width: int
        minimum size of holes to fill in final cell mask
    cell_max_hole_width: int,
        maximum size of holes to fill in final cell mask
    cell_small_obj_w: int
        minimum object size cutoff to remove from final cell mask; likely not required since small objects were removed from cytoplasm mask
    cell_fill_method: str
        method for fill and filter; either "3D" or "slice_by_slice"

    Returns
    -------------
    mask_stack:
        a logical/labels object defining boundaries of nucleus, cell mask, and cytoplasm

    """
    
    cyto_obj = infer_cytoplasm_fromcomposite(in_img,
                                            cyto_weights,
                                            cyto_rescale,
                                            cyto_median_sz,
                                            cyto_gauss_sig,
                                            cyto_mo_method,
                                            cyto_mo_adjust,
                                            cyto_mo_cutoff_size,
                                            cyto_min_hole_w,
                                            cyto_max_hole_w,
                                            cyto_small_obj_w,
                                            cyto_fill_filter_method) 
    
    nuc_obj = infer_nuclei_fromcytoplasm(cyto_obj, 
                                         nuc_min_hole_w,
                                         nuc_max_hole_w,
                                         nuc_fill_method, 
                                         nuc_small_obj_w,
                                         nuc_fill_filter_method)
    
    cell_obj = infer_cellmask_fromcytoplasm(cytoplasm_mask = cyto_obj, 
                                            nucleus_mask = nuc_obj,
                                            min_hole_width = cell_min_hole_width,
                                            max_hole_width = cell_max_hole_width,
                                            small_obj_width = cell_small_obj_width,
                                            fill_filter_method = cell_fill_filter_method)
    
    stack = stack_masks(nuc_mask=nuc_obj, cellmask=cell_obj, cyto_mask=cyto_obj)

    return stack

####################
# Masks B Workflow #
####################
def infer_masks_B(in_img: np.ndarray,
                   cyto_weights: list[int],
                   cyto_rescale: bool,
                   cyto_median_sz: int,
                   cyto_gauss_sig: float,
                   cyto_mo_method: str,
                   cyto_mo_adjust: float,
                   cyto_mo_cutoff_size: int,
                   cyto_min_hole_w: int,
                   cyto_max_hole_w: int,
                   cyto_small_obj_w: int,
                   cyto_fill_filter_method: str,
                   max_nuclei_width: int,
                   nuc_small_obj_width: int,
                   cell_fillhole_max: int,
                   cyto_small_object_width2: int) -> np.ndarray:
    
    """
    Procedure to infer cellmask from linear unmixed input.

    Parameters
    ------------
    in_img: 
        a 3d image containing all the channels
    weights:
        a list of int that corresond to the weights for each channel in the composite; use 0 if a channel should not be included in the composite image
    rescale:
        True = rescale composite
        False = don't rescale composite
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
        determines if fill and filter should be run 'sice-by-slice' or in '3D' 
    cyto_fillhole_max: int
        size of the gaps between the nuclei and cytoplasm (usually small)
    max_nuclei_width: int
        the maximum expected width of the nuclei
    nuc_small_obj_width: int
        width of the any small objects that may be left over after nuclei are selected (i.e., errors in the seg)
    cell_fillhole_max: int
        size of the gaps between the nuclei and cytoplasm (usually small)
    cyto_small_object_width2: int
        size of the gaps between the nuclei and cytoplasm (usually small)

    Returns
    -------------
    mask_stack:
        a three channel np.ndarray constisting of the nucleus, cell and cytoplasm masks (one object per channel)

    """
    cytoplasms = infer_cytoplasm_fromcomposite(in_img, 
                                       cyto_weights,
                                       cyto_rescale,
                                       cyto_median_sz,
                                       cyto_gauss_sig,
                                       cyto_mo_method,
                                       cyto_mo_adjust,
                                       cyto_mo_cutoff_size,
                                       cyto_min_hole_w,
                                       cyto_max_hole_w,
                                       cyto_small_obj_w,
                                       cyto_fill_filter_method)
    
    nuclei_seeds = segment_nuclei_seeds(cytoplasms, 
                              max_nuclei_width, 
                              nuc_small_obj_width)
    
    cellmasks = combine_cytoplasm_and_nuclei(cytoplasms, nuclei_seeds, cell_fillhole_max)
    
    good_CM = select_highest_intensity_cell(in_img, cellmasks, nuclei_seeds)
    
    mask_stack = mask_cytoplasm_nuclei(good_CM, cytoplasms, cyto_small_object_width2)
    
    return mask_stack

####################
# Masks C Workflow #
####################
def infer_masks_C(in_img: np.ndarray,
                   pm_ch: Union[int,None],
                   nuc_ch: Union[int,None],
                   weights_I: list[int],
                   weights_II: list[int],
                   invert_pm_I: bool,
                   invert_pm_II: bool,
                   method_I: str,
                   method_II: str,
                   size_I: int,
                   size_II: int,
                   global_method_I: str,
                   global_method_II: str,
                   cutoff_size_I: int,
                   cutoff_size_II: int,
                   local_adj_I: float,
                   local_adj_II: float,
                   bind_to_pm_I: bool,
                   bind_to_pm_II: bool,
                   thresh_adj_I: float,
                   thresh_adj_II: float,
                   nuc_med_filter_size: int,
                   nuc_gaussian_smoothing_sigma: float,
                   nuc_threshold_factor: float,
                   nuc_thresh_min: float,
                   nuc_thresh_max: float,
                   nuc_hole_min_width: int,
                   nuc_hole_max_width: int,
                   nuc_small_object_width: int,
                   nuc_fill_filter_method: str,
                   nuc_search_img: str,
                   cm_method_I: str,
                   cm_method_II: str,
                   cm_size_I: int,
                   cm_size_II: int,
                   cm_min_hole_width_I: int,
                   cm_min_hole_width_II: int,
                   cm_max_hole_width_I: int,
                   cm_max_hole_width_II: int,
                   cm_small_obj_width_I: int,
                   cm_small_obj_width_II: int,
                   cell_watershed_method: str,
                   cell_min_hole_width: int,
                   cell_max_hole_width: int,
                   cell_method: str,
                   cell_size: int) -> np.ndarray:
    """
    Procedure to infer intermediate masks from linear unmixed input.

    Parameters
    ------------
    in_img: 
        a 3d image containing all the channels
    pm_ch:
        the index of the channel containing your plasma membrane label
    nuc_ch:
        the index of the channel containing your nuclei label
    weights_I:
        a list of int that corresond to the weights for each channel in the first composite; use 0 if a channel should not be included in the composite image
    weights_II:
        a list of int that corresond to the weights for each channel in the second composite; use 0 if a channel should not be included in the composite image
    invert_pm_I:
        True = invert plasma membrane channel in the first composite
        False = do not invert plasma membrane channel in the first composite
    invert_pm_II:
        True = invert plasma membrane channel in the second composite
        False = do not invert plasma membrane channel in the second composite
    method_I:
        which footprint shape to use for the closing algorithm on the first composite. Options include:
        "Ball" (3D closing), "Disk" (2D closing), and "Scharr" (skip closing altogether and apply Scharr edge detection).
    method_II:
        which footprint shape to use for the closing algorithm on the second composite. Options include:
        "Ball" (3D closing), "Disk" (2D closing), and "Scharr" (skip closing altogether and apply Scharr edge detection).
    size_I:
        size of the footprint used in closing algorithm on the first composite, this value is disregarded if method_I == "Scharr"
    size_II:
        size of the footprint used in closing algorithm on the second composite, this value is disregarded if method_II == "Scharr"
    global_method_I:
         which method to use for calculating global threshold applied to composite_mask_I (MO). Options include:
         "triangle" (or "tri"), "median" (or "med"), and "ave_tri_med" (or "ave").
         "ave" refers the average of "triangle" threshold and "mean" threshold.
    global_method_II:
         which method to use for calculating global threshold applied to composite_mask_II (MO). Options include:
         "triangle" (or "tri"), "median" (or "med"), and "ave_tri_med" (or "ave").
         "ave" refers the average of "triangle" threshold and "mean" threshold.
    cutoff_size_I: 
        Masked Object threshold `size_min`; minimum size of of object to advanced to the local Otsu thresholding
        step in the Masked Object thresholding of composite_mask_I.
    cutoff_size_II: 
        Masked Object threshold `size_min`; minimum size of of object to advanced to the local Otsu thresholding
        step in the Masked Object thresholding of composite_mask_II.
    local_adj_I: 
        Masked Object threshold `local_adjust`, proportion applied to the local Otsu threshold (composite_mask_I MO thresholding)
    local_adj_II: 
        Masked Object threshold `local_adjust`, proportion applied to the local Otsu threshold (composite_mask_II MO thresholding)
    bind_to_pm_I:
        True = restrict the resulting mask_I_segmentation to the thresholded plasma membrane
        False = do not restrict the resulting mask_I_segmentation to the thresholded plasma membrane
    bind_to_pm_II:
        True = restrict the resulting mask_II_segmentation to the thresholded plasma membrane
        False = do not restrict the resulting mask_II_segmentation to the thresholded plasma membrane
    thresh_adj_I:
        the proportion applied to the Otsu threshold of the plasma membrane (disregarded if bind_to_pm_I == False)
    thresh_adj_II:
        the proportion applied to the Otsu threshold of the plasma membrane (disregarded if bind_to_pm_II == False)
    nuc_med_filter_size: 
        width of median filter for nuclei signal
    nuc_gaussian_smoothing_sigma: 
        sigma for gaussian smoothing of nuclei signal
    nuc_threshold_factor:
        adjustment to make to the intial local threshold of the nucleus singal;
        intital threshold is derived from Li's minimum cross entropy method
    nuc_thresh_min:
        minimum bound of the local threshold of the nucleus
    nuc_thresh_max:
        maximum bound of the local threshold of the nucleus
    nuc_hole_min_width: 
        the minimum hole width to be filled in the nucleus post-thresholding
    nuc_hole_max_width:
        the maximum hole width to be filled in the nucleus post-thresholding
    nuc_small_object_width:
        minimum object size cutoff for nucleus post-thresholding
    nuc_fill_filter_method:
        determines if fill and filter should be run 'sice-by-slice' or in '3D'
    nuc_search_img:
        the segmentation used to select the nucleus based on greatest overlap. Options include:
        "Img 5" (mask_I_segmentation) and "Img 6" (mask_II_segmentation)
    cm_method_I:
        which footprint shape to use for the nucleus dilation before combination with the mask_I_segmentation. Options include:
        "Ball" (3D dilation), "Disk" (2D dilation), and "None" (skip dilation altogether).
    cm_method_II:
        which footprint shape to use for the nucleus dilation before combination with the mask_II_segmentation. Options include:
        "Ball" (3D dilation), "Disk" (2D dilation), and "None" (skip dilation altogether).
    cm_size_I:
        size of the footprint used in dilation of the nucleus object, this value is disregarded if cm_method_I == "None"
    cm_size_II:
        size of the footprint used in dilation of the nucleus object, this value is disregarded if cm_method_II == "None"
    cm_min_hole_width_I: 
        the minimum hole width to be filled in the cellmask_I_segmentation
    cm_min_hole_width_II: 
        the minimum hole width to be filled in the cellmask_II_segmentation
    cm_max_hole_width_I: 
        the maximum hole width to be filled in the cellmask_I_segmentation
    cm_max_hole_width_II: 
        the maximum hole width to be filled in the cellmask_II_segmentation
    cm_small_obj_width_I:
        minimum object size cutoff for the cellmask_I_segmentation
    cm_small_obj_width_II:
        minimum object size cutoff for the cellmask_II_segmentation
    cell_watershed_method:
        determines if the watershed should be run 'sice-by-slice' or in '3D'
    cell_min_hole_width: 
        the minimum hole width to be filled in the cellmask object
    cell_max_hole_width: 
        the maximum hole width to be filled in the cellmask object
    cell_method:
        which footprint shape to use for the cellmask closing (dilation -> erosion). Options include:
        "Ball" (3D closing), "Disk" (2D closing), and "None" (skip closing altogether).
    cell_size:
        size of the footprint used in closing of the cellmask object, this value is disregarded if cm_method_I == "None"
    
    Returns
    -------------
    mask_stack:
        a two channel np.ndarray constisting of the nucleus and cell (one object per channel)

    """

    ##########################
    # get intermediate masks #
    ##########################

    mask_I_segmentation, mask_II_segmentation = infer_intermediate_masks(in_img,
                           pm_ch,
                           weights_I,
                           weights_II,
                           invert_pm_I,
                           invert_pm_II,
                           method_I,
                           method_II,
                           size_I,
                           size_II,
                           global_method_I,
                           global_method_II,
                           cutoff_size_I,
                           cutoff_size_II,
                           local_adj_I,
                           local_adj_II,
                           bind_to_pm_I,
                           bind_to_pm_II,
                           thresh_adj_I,
                           thresh_adj_II)
    
    #################
    # infer nucleus #
    #################
    
    nuc_obj = infer_nucleus_masks_C(in_img,
                           nuc_ch,
                           mask_I_segmentation,
                           mask_II_segmentation,
                           nuc_med_filter_size,
                           nuc_gaussian_smoothing_sigma,
                           nuc_threshold_factor,
                           nuc_thresh_min,
                           nuc_thresh_max,
                           nuc_hole_min_width,
                           nuc_hole_max_width,
                           nuc_small_object_width,
                           nuc_fill_filter_method,
                           nuc_search_img)
    
    ##################
    # infer cellmask #
    ##################
    
    cellmask_obj = infer_cellmask_masks_C(nuc_obj,
                                    mask_I_segmentation,
                                    mask_II_segmentation,
                                    cm_method_I,
                                    cm_method_II,
                                    cm_size_I,
                                    cm_size_II,
                                    cm_min_hole_width_I,
                                    cm_min_hole_width_II,
                                    cm_max_hole_width_I,
                                    cm_max_hole_width_II,
                                    cm_small_obj_width_I,
                                    cm_small_obj_width_II,
                                    cell_watershed_method,
                                    cell_min_hole_width,
                                    cell_max_hole_width,
                                    cell_method,
                                    cell_size)
    
    ###################
    ### stack masks ###
    ###################
    mask_stack = stack_masks(nuc_mask=nuc_obj, 
                            cellmask=cellmask_obj)
    
    return mask_stack

def infer_masks_D(in_img: np.ndarray,
                   pm_ch: Union[int,None],
                   nuc_ch: Union[int,None],
                   weights: list[int],
                   median_sz: int, 
                   gauss_sig: float,
                   thresh_factor: float,
                   thresh_min: float,
                   thresh_max: float,
                   min_hole_w: int,
                   max_hole_w: int,
                   small_obj_w: int,
                   fill_filter_method: str,
                   invert_pm: bool,
                   cell_method: str,
                   hole_min: int,
                   hole_max: int,
                   fill_2d: bool):
    
    """
    Procedure to infer intermediate masks from linear unmixed input.

    Parameters
    ------------
    in_img: 
        a 3d image containing all the channels
    pm_ch:
        the index of the channel containing your plasma membrane label
    nuc_ch:
        the index of the channel containing your nuclei label
    weights:
        a list of int that corresond to the weights for each channel in the composite; use 0 if a channel should not be included in the composite image
    median_sz: int
        width of median filter for nuclei signal
    gauss_sig: float
        sigma for gaussian smoothing of nuclei signal
    thresh_factor: float
        adjustment factor for log Li nuclei threholding
    thresh_min: float
        abs min threhold for log Li nuclei threholding
    thresh_max: float
        abs max threhold for log Li nuclei threholding
    min_hole_w: int
        hole filling minimum for nuclei post-processing
    max_hole_w: int
        hole filling cutoff for nuclei post-processing
    small_obj_w: int
        minimum object size cutoff for nuclei post-processing
    fill_filter_method:
        determines if fill and filter should be run 'sice-by-slice' or in '3D' (nucleus)
    invert_pm:
        True = invert plasma membrane channel in the composite
        False = do not invert plasma membrane channel in the composite
    cell_method:
        The type of watershedding method to perform for the cellmask. Options include:
        "3D" (3D watershedding), "slice-by-slice" (2D watershedding).
    hole_min: 
        the minimum hole width to be filled in the cellmask
    hole_max: 
        the maximum hole width to be filled in the cellmask
    fill_2d:
        True = fill cellmask holes in 2D (slice by slice)
        False = fill cellmask holes in 3D
    
    Returns
    -------------
    mask_stack:
        a two channel np.ndarray constisting of the nucleus and cell (one object per channel)

    """

    ##########################
    # get nuclei labels #
    ##########################

    nuclei_labels = infer_nuclei_fromlabel(in_img,
                                nuc_ch,
                                median_sz,
                                gauss_sig,
                                thresh_factor,
                                thresh_min,
                                thresh_max,
                                min_hole_w,
                                max_hole_w,
                                small_obj_w,
                                fill_filter_method)
    
    ##########################
    # get cellmask #
    ##########################

    cellmask_obj = infer_cellmask_masks_D(in_img,
                                    pm_ch,
                                    weights,
                                    invert_pm,
                                    nuclei_labels,
                                    cell_method,
                                    hole_min,
                                    hole_max,
                                    fill_2d)
    
    ###################
    ### stack layers ##
    ###################

    stack = stack_layers(nuclei_labels = nuclei_labels,
                         cellmask = cellmask_obj)
    
    return stack