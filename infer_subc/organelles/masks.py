from typing import Union

import numpy as np

from infer_subc.core.img import *
from infer_subc.organelles.cellmask import infer_cellmask_fromcomposite
from infer_subc.organelles.nuclei import infer_nuclei_fromlabel
from infer_subc.organelles.cytoplasm import infer_cytoplasm


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
        minimu object size cutoff for nuclei post-processing for nuclei channel
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
        minimu object size cutoff for cellmask signal post-processing
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