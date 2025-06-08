import numpy as np
from skimage.morphology import disk, ball, binary_closing, closing, binary_dilation, dilation, binary_erosion
from infer_subc.core.img import fill_and_filter_linear_size, get_interior_labels, get_max_label, inverse_log_transform, log_transform, masked_inverted_watershed, masked_object_thresh, min_max_intensity_normalization, select_channel_from_raw, threshold_otsu_log
from aicssegmentation.core.utils import hole_filling
from infer_subc.organelles.nuclei import infer_nuclei_fromlabel
from infer_subc.organelles.cellmask import non_linear_cellmask_transform
from infer_subc.core.img import label_bool_as_uint16
from typing import Union

def membrane_composite(in_img: np.ndarray,
                       weight_ch0: int = 0,
                       weight_ch1: int = 0,
                       weight_ch2: int = 0,
                       weight_ch3: int = 0,
                       weight_ch4: int = 0,
                       weight_ch5: int = 0,
                       weight_ch6: int = 0,
                       weight_ch7: int = 0,
                       weight_ch8: int = 0,
                       weight_ch9: int = 0,
                       Invert_PM: bool = False,
                       PM_Channel: int = 0,
                     ):
    weights = [weight_ch0,
               weight_ch1,
               weight_ch2,
               weight_ch3,
               weight_ch4,
               weight_ch5,
               weight_ch6,
               weight_ch7,
               weight_ch8,
               weight_ch9]
    raw_img = in_img.copy()
    out_img = np.zeros_like(in_img[0]).astype(np.double)
    if Invert_PM:
        raw_img[PM_Channel] = np.max(in_img[PM_Channel])-in_img[PM_Channel]
    for ch, weight in enumerate(weights):
        if weight > 0:
            out_img+=weight*min_max_intensity_normalization(raw_img[ch])
    return out_img

def masked_object_thresh_bind_pm(raw_img: np.ndarray,
                                 in_img: np.ndarray,
                                 Global_Method: str,
                                 Cutoff_Size: int,
                                 Local_Adjust: float,
                                 Bind_to_PM: bool,            
                                 PM_Channel: int,
                                 Thresh_Adj: float):
    threshed = masked_object_thresh(in_img, Global_Method, Cutoff_Size, Local_Adjust)

    if Bind_to_PM:
        pm_img = select_channel_from_raw(raw_img, PM_Channel)
        pm_image, d = log_transform(pm_img.copy())
        pm_thresh = threshold_otsu_log(pm_image) 
        invert_pm_obj = np.invert(pm_img > (inverse_log_transform(pm_thresh, d) * Thresh_Adj))
        mask = np.zeros_like(invert_pm_obj)
        mask[(invert_pm_obj == threshed) & 
                 (invert_pm_obj == True) & 
                 (threshed == True)] = 1
        out_img = closing(mask)
    else:
        out_img = threshed.copy()

    return out_img

def close_and_filter(in_img: np.ndarray, 
          Method: str, 
          Size: int):
    out_img = np.zeros_like(in_img)
    if Method == "Ball":
        if Size > 0:
            fp = ball(Size)
            out_img = closing(in_img.copy(), footprint=fp)
        else:
            out_img = closing(in_img.copy())
    elif Method == "Disk":
        if Size > 0:
            fp = disk(Size)
            for z in range(len(in_img)):
                out_img[z] = closing(in_img[z].copy(), footprint=fp)
        else:
            for z in range(len(in_img)):
                out_img[z] = closing(in_img[z].copy())
    elif Method == "Scharr":
        out_img = non_linear_cellmask_transform(in_img)
    elif Method == "None":
        out_img = in_img
    return out_img

def find_nuc(raw_img: np.ndarray, 
             in_img_A: np.ndarray,
             in_img_B: np.ndarray,
             Nuc_Channel: int,
             Median_Size: int,
             Gauss_Sigma: float,
             Thresh_Factor: float,
             Thresh_Min: float,
             Thresh_Max: float,
             Min_Hole_Width: int,
             Max_Hole_Width: int,
             Small_Obj_Width: int,
             Fill_Filter_Method: str,
             Search_Img: str
             ):
    nuc_img = infer_nuclei_fromlabel(in_img=raw_img, 
                                     nuc_ch=Nuc_Channel,
                                     median_sz=Median_Size, 
                                     gauss_sig=Gauss_Sigma,
                                     thresh_factor=Thresh_Factor,
                                     thresh_min=Thresh_Min,
                                     thresh_max=Thresh_Max,
                                     min_hole_w=Min_Hole_Width,
                                     max_hole_w=Max_Hole_Width,
                                     small_obj_w=Small_Obj_Width,
                                     fill_filter_method=Fill_Filter_Method)
    if Search_Img == "Img 5":
        in_img = in_img_A
    elif Search_Img == "Img 6":
        in_img = in_img_B
    keep_nuc = get_max_label((in_img), dilation(nuc_img))
    out_img = np.zeros_like(nuc_img)
    out_img[nuc_img == keep_nuc] = 1
    return out_img

def mix_nuc_and_fill(nuc: np.ndarray, 
            in_img: np.ndarray, 
            Method: str, 
            Size: int,
            Min_Hole_Width: int,
            Max_Hole_Width: int,
            Small_Obj_Width: int):
    single_nuc = nuc
    out_img = in_img.copy()
    if Method == "Disk":
        if Size > 0:    
            nuc_fp = disk(Size)
            for z in range(len(in_img)):
                out_img[z] += binary_dilation(single_nuc.astype(bool)[z], footprint=nuc_fp)
        else:
            for z in range(len(in_img)):
                out_img[z] += binary_dilation(single_nuc.astype(bool)[z])
    elif Method == "Ball":
        if Size > 0:  
            nuc_fp = ball(Size)
            out_img += binary_dilation(single_nuc.astype(bool), footprint=nuc_fp)
        else:
            out_img += binary_dilation(single_nuc.astype(bool))
    elif Method == "None":
        out_img += single_nuc.astype(bool)
    filled = fill_and_filter_linear_size(out_img,
                                         hole_min=Min_Hole_Width,
                                         hole_max=Max_Hole_Width,
                                         min_size=Small_Obj_Width)
    return filled

def fill_and_bind(raw_img: np.ndarray, 
                  in_img: np.ndarray, 
                  Min_Hole_Width: int, 
                  Max_Hole_Width: int,
                  Small_Obj_Width: int,
                  PM_Channel: int,
                  Thresh_Adj: float,
                  Bind_to_PM: bool):
    filled = fill_and_filter_linear_size(in_img,
                                         hole_min=Min_Hole_Width,
                                         hole_max=Max_Hole_Width,
                                         min_size=Small_Obj_Width)
    if Bind_to_PM:
        pm_img = select_channel_from_raw(raw_img, PM_Channel)
        pm_image, d = log_transform(pm_img.copy())
        pm_thresh = threshold_otsu_log(pm_image) 
        invert_pm_obj = np.invert(pm_img > (inverse_log_transform(pm_thresh, d) * Thresh_Adj))
        mask = np.zeros_like(invert_pm_obj)
        mask[(binary_dilation(invert_pm_obj) > 0) & 
             (filled > 0)] = 1
        out_img = binary_closing(mask)
    else:
        out_img = filled.copy()
    return out_img

def close_and_fill(in_img: np.ndarray, 
                   Min_Hole_Width: int,
                   Max_Hole_Width: int,
                   Method: str,
                   Size: int):
    if Method == "Disk":
        if Size > 0:
            close_fp = disk(Size)
            for z in range(len(in_img)):
                in_img[z] = binary_dilation(in_img.copy()[z], footprint=close_fp)
        else:
            for z in range(len(in_img)):
                in_img[z] = binary_dilation(in_img.copy()[z])
    elif Method == "Ball":
        if Size > 0:
            close_fp = ball(Size)
            in_img = binary_dilation(in_img.copy(), footprint=close_fp)
        else:
            in_img = binary_dilation(in_img.copy())
    in_img = hole_filling(in_img.copy(), hole_min=Min_Hole_Width, hole_max=Max_Hole_Width, fill_2d=True)
    if Method == "Disk":
        if Size > 0:
            for z in range(len(in_img)):
                in_img[z] = binary_erosion(in_img.copy()[z], footprint=close_fp)
        else:
            for z in range(len(in_img)):
                in_img[z] = binary_erosion(in_img.copy()[z])
    elif Method == "Ball":
        if Size > 0:
            in_img = binary_erosion(in_img.copy(), footprint=close_fp)
        else:
            in_img = binary_erosion(in_img.copy())
    return in_img

# def double_watershed(nuc: np.ndarray,
#                      raw_img_A: np.ndarray,
#                      raw_img_B: np.ndarray,
#                      thresh_img_A: np.ndarray,
#                      thresh_img_B: np.ndarray,
#                      Watershed_Method: str,
#                      Min_Hole_Width: int,
#                      Max_Hole_Width: int,
#                      Method: str,
#                      Size: int):
    
#     choose_nuc = nuc

#     cm_A = masked_inverted_watershed(raw_img_A, 
#                                      choose_nuc, 
#                                      thresh_img_A,
#                                      method=Watershed_Method)
#     cm_B = masked_inverted_watershed(raw_img_B, 
#                                      choose_nuc, 
#                                      thresh_img_B,
#                                      method=Watershed_Method)
    
#     cm_combo = cm_A.astype(bool) + cm_B.astype(bool)
#     cm_out = close_and_fill(in_img=cm_combo,
#                              Min_Hole_Width=Min_Hole_Width,
#                              Max_Hole_Width=Max_Hole_Width,
#                              Method=Method,
#                              Size=Size)
    
#     out_img = np.stack([nuc, cm_out], axis=0)

#     return out_img

def double_watershed(nuc: np.ndarray,
                     raw_img_A: np.ndarray,
                     raw_img_B: np.ndarray,
                     thresh_img_A: np.ndarray,
                     thresh_img_B: np.ndarray,
                     Watershed_Method: str,
                     Min_Hole_Width: int,
                     Max_Hole_Width: int,
                     Method: str,
                     Size: int):
    
    choose_nuc = nuc

    cm_A = masked_inverted_watershed(raw_img_A, 
                                     choose_nuc, 
                                     thresh_img_A,
                                     method=Watershed_Method)
    cm_B = masked_inverted_watershed(raw_img_B, 
                                     choose_nuc, 
                                     thresh_img_B,
                                     method=Watershed_Method)
    
    cm_combo = cm_A.astype(bool) + cm_B.astype(bool)
    cm_out = close_and_fill(in_img=cm_combo,
                             Min_Hole_Width=Min_Hole_Width,
                             Max_Hole_Width=Max_Hole_Width,
                             Method=Method,
                             Size=Size)

    return cm_out

def invert_pm_watershed(raw_img: np.ndarray, nuc_labels: np.ndarray, PM_Channel: int, Method: str):
    pm_img = select_channel_from_raw(raw_img, PM_Channel)
    invert_pm_img = abs(np.max(pm_img) - pm_img)
    out_img = masked_inverted_watershed(img_in=min_max_intensity_normalization(invert_pm_img),
                                     markers=nuc_labels,
                                     mask=None,
                                     method=Method)
    return out_img

def choose_cell(composite_img: np.ndarray, nuc_labels: np.ndarray, watershed: np.ndarray):
    target_labels = get_interior_labels(nuc_labels)
    keep_label = get_max_label(composite_img, 
                           watershed,
                           target_labels=target_labels)
    cellmask_out = np.zeros_like(watershed)
    cellmask_out[watershed == keep_label] = 1
    return cellmask_out

##########################
# infer_intermediate_masks
# alternative workflow "c"
##########################

def infer_intermediate_masks(in_img: np.ndarray,
                           pm_ch: Union[int,None],
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
                           thresh_adj_II: float):
 
    """
    Procedure to infer intermediate masks from linear unmixed input.

    Parameters
    ------------
    in_img: 
        a 3d image containing all the channels
    pm_ch:
        the index of the channel containing your plasma membrane label
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
    
    

    Returns
    -------------
    mask_I and mask_II segmentation:
        two logical/labels objects that will be used to find the cellmask

    """

    ###################
    # EXTRACT
    ###################

    struct_img_raw_I = membrane_composite(in_img, *weights_I, invert_pm_I, pm_ch)

    struct_img_raw_II = membrane_composite(in_img, *weights_II, invert_pm_II, pm_ch)

    ###################
    # PRE_PROCESSING
    ###################

    composite_mask_I = close_and_filter(struct_img_raw_I, method_I, size_I)

    composite_mask_II = close_and_filter(struct_img_raw_II, method_II, size_II)

    ###################
    # CORE_PROCESSING
    ###################

    mask_I_segmentation = masked_object_thresh_bind_pm(in_img,
                                                       composite_mask_I,
                                                       global_method_I,
                                                       cutoff_size_I,
                                                       local_adj_I,
                                                       bind_to_pm_I,
                                                       pm_ch,
                                                       thresh_adj_I)
    
    mask_II_segmentation = masked_object_thresh_bind_pm(in_img,
                                                       composite_mask_II,
                                                       global_method_II,
                                                       cutoff_size_II,
                                                       local_adj_II,
                                                       bind_to_pm_II,
                                                       pm_ch,
                                                       thresh_adj_II)
    
    return label_bool_as_uint16(mask_I_segmentation),label_bool_as_uint16(mask_II_segmentation)
