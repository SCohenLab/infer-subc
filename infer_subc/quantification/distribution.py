import itertools
from pathlib import Path
import time

import numpy as np
import pandas as pd

from skimage.measure import regionprops_table, regionprops, mesh_surface_area, marching_cubes, label

import centrosome.cpmorphology
import centrosome.propagate
import centrosome.zernike

from infer_subc.core.img import apply_mask
from infer_subc.utils.batch import list_image_files, find_segmentation_tiff_files
from infer_subc.core.file_io import read_czi_image, read_tiff_image
from infer_subc.quantification.batch import append_atomic_csv, load_existing_keys_csv
from typing import Tuple, Any, Union, List

# from scipy.ndimage import maximum_position, center_of_mass
from scipy.ndimage import sum as ndi_sum
from scipy.sparse import coo_matrix


###################################
### XY DISTRIBUTION
###################################

### USED ###
def create_masked_sum_projection(img_in:np.ndarray, mask:Union[np.ndarray, None]=None, to_bool:bool=True) -> np.ndarray:
    """
    Parameters:
    ----------
    img_in:
        3D (ZYX) np.ndarray that will be summed along the Z axis
    mask:
        Optional - mask of the region you want to include in the final sum projection
    to_bool:
        True = input image is created in a boolean image before sum projection (useful for segmentation images where each object is coded as a unique ID number; like after skimage.segmentation.label())
        False = original input image is used for the sum projection
    """
    img_out = img_in.astype(bool) if to_bool else img_in
    if mask is not None:
        img_out = apply_mask(img_out, mask)
    
    return img_out.sum(axis=0)

def get_normalized_distance_and_mask(labels: np.ndarray, 
                                      center_objects: Union[np.ndarray, None], 
                                      center_on: bool,
                                      intres: Union[bool, None] = False):
    """
    helper for radial distribution
    Parameters:
    ----------
    labels:
        2D (YX) np.ndarray - normally the result of a binary ZYX segmentation of the cell mask after a sum projection across the Z dimension.
        If labels does not contain a true mask/background (e.g., labels.min() != 0), the distance to the outer edge is not used to
        define the normalized_distance values, although d_to_edge may still be computed (in particular when center_objects is None).
    center_object:
        2D (YX) np.ndarray - normally the result of a binary ZYX segmentation of the nucleus after a sum projection across the Z dimension.
        If no centering object is included, the center of the labels will be used.
    center_on:
        True = the center of the centering object will be used as the starting point to calculate the distance from the center
        False = the edge of the centering object will be used as the starting point to calculate the distance from the center
    intres:
        True = output 4 default objects in addition to d_to_edge and d_from_center np.ndarrays
        False = only output 4 default objects
    
    Output:
    ----------
    normalized_distance:
        2D (YX) np.ndarray with intensity values representing the distance between the edge of the "labels" and the centering object.
        More specifically for the centering object, either the edge or the centermost point is used, depending on the center_on
        parameter. If there is no centering object, the values will represent the distance from the edge of the "labels" object. 
    good_mask:
        mask of the areas that were included in the normalized_distance output
    i_center: If center_objects *is not* None: i (Y) coordinate of the centermost point of the centering object
              If center_objects *is* None: the i (Y) coordinate of the innermost (distance from the edge) point of the "labels" input
    j_center: If center_objects *is not* None: j (X) coordinate of the centermost point of the centering object
              If center_objects *is* None: the j (X) coordinate of the innermost (distance from the edge) point of the "labels" input
    d_to_edge:
        2D (YX) np.ndarray with intensity values representing the distance from the edge of the "labels" object (if labels is not None).
    d_from_center:
        2D (YX) np.ndarray with intensity values representing the distance from the centermost point of the centering object,
        or the edge of the centering object (depending on the value of center_on).

    """

    # apply a euclidian distance transform for the cellmask projection (if one exists); brightness represents the distance from the edge of the cell
    # First case, there exists a "true" mask object (labels is expected to be BINARY)
    if labels.min() == 0:
        true_mask = True
        d_to_edge = centrosome.cpmorphology.distance_to_edge(labels)
    # Secondary case, no "true" mask object
    else:
        true_mask = False 
        print('labels/mask input represents entire image')
        if center_objects is None:
            # pad the labels array by zeros so that the edge pixels is detected as boundary pixels
            padded_labels = np.pad(labels, pad_width = 1)
            # the cropped image is then set to d_to_edge (only used to find centermost point of image frame)
            d_to_edge = centrosome.cpmorphology.distance_to_edge(padded_labels)[1:-1,1:-1]
        else:
            d_to_edge = None
            
            

    if center_objects is not None:
        
        # this lists the pixel counts for each cell mask in the image based on the number of unique centering objects
        center_labels = label(center_objects)
        pixel_counts = centrosome.cpmorphology.fixup_scipy_ndimage_result(ndi_sum(np.ones(center_labels.shape), 
                                                                                  center_labels, 
                                                                                  np.arange(1, np.max(center_labels) + 1, dtype=np.int32)))
        good = pixel_counts > 0
        i, j = (centrosome.cpmorphology.centers_of_labels(center_labels) + 0.5).astype(int)
        ig = i[good]
        jg = j[good]
        lg = np.arange(1, len(i) + 1)[good]
        
        if center_on:  # Reduce the propagation labels to the centers of the centering objects
            center_labels = np.zeros(center_labels.shape, int)
            center_labels[ig, jg] = lg

        cl, d_from_center = centrosome.propagate.propagate(np.zeros(center_labels.shape), center_labels, labels != 0, 1)
        cl[labels == 0] = 0

        missing_mask = (labels != 0) & (cl == 0)
        missing_labels = np.unique(labels[missing_mask])
        
        if len(missing_labels):
            print("how did we have missing labels?")
            all_centers = centrosome.cpmorphology.centers_of_labels(labels)
            missing_i_centers, missing_j_centers = all_centers[:, missing_labels-1]
            di = missing_i_centers[:, np.newaxis] - ig[np.newaxis, :]
            dj = missing_j_centers[:, np.newaxis] - jg[np.newaxis, :]
            missing_best = lg[np.argsort(di * di + dj * dj)[:, 0]]
            best = np.zeros(np.max(labels) + 1, int)
            best[missing_labels] = missing_best
            cl[missing_mask] = best[labels[missing_mask]]

            iii, jjj = np.mgrid[0 : labels.shape[0], 0 : labels.shape[1]]
            di = iii[missing_mask] - i[cl[missing_mask] - 1]
            dj = jjj[missing_mask] - j[cl[missing_mask] - 1]
            d_from_center[missing_mask] = np.sqrt(di * di + dj * dj)

        good_mask = cl > 0
            
    else:
        # i, j = centrosome.cpmorphology.maximum_position_of_labels(d_to_edge, labels, [1])
        i, j = centrosome.cpmorphology.maximum_position_of_labels(d_to_edge, labels, [1])
        # delete d_to_edge if no true mask exists
        if not true_mask:
            d_to_edge = None 
        center_labels = np.zeros(labels.shape, int)
        center_labels[i, j] = labels[i, j]
        colors = centrosome.cpmorphology.color_labels(labels)
        ncolors = np.max(colors)
        d_from_center = np.zeros(labels.shape)
        cl = np.zeros(labels.shape, int)

        for color in range(1, ncolors + 1):
            mask = colors == color
            # There is no Z height if we literally have flattened the image
            l, d = centrosome.propagate.propagate( np.zeros(center_labels.shape), center_labels, mask, 1)
            d_from_center[mask] = d[mask]
            cl[mask] = l[mask]

        good_mask = cl > 0

    # creating an object equal to the cellmask_proj with all pixel values equal to the Y coordinate value (here called 'i') or X coordinate (here called 'j')
    # then creating normalized distance out from center to edge of cell (if mask exists)
    i_center = np.zeros(cl.shape)
    i_center[good_mask] = i[cl[good_mask] - 1]

    j_center = np.zeros(cl.shape)
    j_center[good_mask] = j[cl[good_mask] - 1]

    normalized_distance = np.zeros(labels.shape)
    total_distance = d_from_center + d_to_edge if true_mask else d_from_center

    # Normalize the total distance
    if true_mask:
        normalized_distance[good_mask] = d_from_center[good_mask] / (total_distance[good_mask] + 0.001)
    else:
        normalized_distance[good_mask] = d_from_center[good_mask] / (d_from_center.max() + 0.001)
    
    # include d_to_edge and d_from_center?
    if intres:
        return normalized_distance, good_mask, i_center, j_center, d_to_edge, d_from_center
    else:
        return normalized_distance, good_mask, i_center, j_center
# def get_normalized_distance_and_mask(labels: np.ndarray, 
#                                       center_objects: Union[np.ndarray, None], 
#                                       center_on: bool,
#                                       intres: Union[bool, None] = False):
#     """
#     helper for radial distribution
#     Parameters:
#     ----------
#     labels:
#         2D (YX) np.ndarray - normally the result of a binary ZYX segmentation of the cell mask after a sum projection across the Z dimension
#     center_object:
#         2D (YX) np.ndarray - normally the result of a binary ZYX segmentation of the nucleus after a sum projection across the Z dimension.
#         If no centering object is included, the center of the labels will be used.
#     center_on:
#         True = the center of the centering object will be used as the starting point to calculate the distance from the center
#         False = the edge of the centering object will be used as the starting point to calculate the distance from the center
#     intres:
#         True = output 4 default objects as wells as d_to_edge and d_from_center np.ndarrays
#         False = output 4 default objects
    
#     Output:
#     ----------
#     normalized_distance:
#         2D (YX) np.ndarray with intensity values representing the distance between the edge of the "labels" and the centering object.
#         More specifically for the centering object, either the edge or the centermost point is used, depending on the center_on
#         parameter. If there is no centering object, the values will represent the distance from the edge of the "labels" object.
#     good_mask:
#         mask of the areas that were included in the normalized_distance output
#     i_center: If center_objects *is not* None: i coordinate of the centermost point of the centering object
#               If center_objects *is* None: the i coordinate of the innermost (distance from the edge) point of the "labels" input
#     j_center: If center_objects *is not* None: j coordinate of the centermost point of the centering object
#               If center_objects *is* None: the j coordinate of the innermost (distance from the edge) point of the "labels" input
#     d_to_edge:
#         2D (YX) np.ndarray with intensity values representing the distance from the edge of the "labels" object.
#     d_from_center:
#         2D (YX) np.ndarray with intensity values representing the distance from the centermost point of the centering object,
#         or the edge of the centering object (depending on the value of center_on).

#     """

#     d_to_edge = centrosome.cpmorphology.distance_to_edge(labels)

#     if center_objects is not None:
#         center_labels = label(center_objects)
#         pixel_counts = centrosome.cpmorphology.fixup_scipy_ndimage_result(ndi_sum(np.ones(center_labels.shape), 
#                                                                                   center_labels, 
#                                                                                   np.arange(1, np.max(center_labels) + 1, dtype=np.int32)))
#         good = pixel_counts > 0
#         i, j = (centrosome.cpmorphology.centers_of_labels(center_labels) + 0.5).astype(int)
#         ig = i[good]
#         jg = j[good]
#         lg = np.arange(1, len(i) + 1)[good]
        
#         if center_on:  # Reduce the propagation labels to the centers of the centering objects
#             center_labels = np.zeros(center_labels.shape, int)
#             center_labels[ig, jg] = lg

#         cl, d_from_center = centrosome.propagate.propagate(np.zeros(center_labels.shape), center_labels, labels != 0, 1)
#         cl[labels == 0] = 0

#         missing_mask = (labels != 0) & (cl == 0)
#         missing_labels = np.unique(labels[missing_mask])
        
#         if len(missing_labels):
#             print("how did we have missing labels?")
#             all_centers = centrosome.cpmorphology.centers_of_labels(labels)
#             missing_i_centers, missing_j_centers = all_centers[:, missing_labels-1]
#             di = missing_i_centers[:, np.newaxis] - ig[np.newaxis, :]
#             dj = missing_j_centers[:, np.newaxis] - jg[np.newaxis, :]
#             missing_best = lg[np.argsort(di * di + dj * dj)[:, 0]]
#             best = np.zeros(np.max(labels) + 1, int)
#             best[missing_labels] = missing_best
#             cl[missing_mask] = best[labels[missing_mask]]

#             iii, jjj = np.mgrid[0 : labels.shape[0], 0 : labels.shape[1]]
#             di = iii[missing_mask] - i[cl[missing_mask] - 1]
#             dj = jjj[missing_mask] - j[cl[missing_mask] - 1]
#             d_from_center[missing_mask] = np.sqrt(di * di + dj * dj)

#         good_mask = cl > 0
            
#     else:
#         # i, j = centrosome.cpmorphology.maximum_position_of_labels(d_to_edge, labels, [1])
#         i, j = centrosome.cpmorphology.maximum_position_of_labels(d_to_edge, labels, [1])
#         center_labels = np.zeros(labels.shape, int)
#         center_labels[i, j] = labels[i, j]
#         colors = centrosome.cpmorphology.color_labels(labels)
#         ncolors = np.max(colors)
#         d_from_center = np.zeros(labels.shape)
#         cl = np.zeros(labels.shape, int)

#         for color in range(1, ncolors + 1):
#             mask = colors == color
#             # There is no Z height if we literally have flattened the image
#             l, d = centrosome.propagate.propagate( np.zeros(center_labels.shape), center_labels, mask, 1)
#             d_from_center[mask] = d[mask]
#             cl[mask] = l[mask]

#         good_mask = cl > 0

#     i_center = np.zeros(cl.shape)
#     i_center[good_mask] = i[cl[good_mask] - 1]

#     j_center = np.zeros(cl.shape)
#     j_center[good_mask] = j[cl[good_mask] - 1]

#     normalized_distance = np.zeros(labels.shape)
#     total_distance = d_from_center + d_to_edge
#     normalized_distance[good_mask] = d_from_center[good_mask] / (total_distance[good_mask] + 0.001)
    
#     # include d_to_edge and d_from_center?
#     if intres:
#         return normalized_distance, good_mask, i_center, j_center, d_to_edge, d_from_center
#     else:
#         return normalized_distance, good_mask, i_center, j_center
    

### USED ###
def get_concentric_distribution(
        mask_proj: np.ndarray,
        mask_name: str,
        centering_proj: np.ndarray,
        obj_proj: np.ndarray,
        obj_name: str,
        bin_count: int,
        center_on: bool = False,
        keep_center_as_bin: bool = True,
        scale: Union[tuple, None]=None):
    """
    Based on CellProfiler's measureobjectintensitydistribution. Measure the distribution of segmented objects within a masked area. 
    In our case, we will usually utilize this function to measure the amount of an organelle within the cell.
    Radial bins are created out from a center point, usually the nucleus edge.

    
    Parameters
    ------------
    mask_proj: np.ndarray
        a sum projection of the region you want to measure the distribution from where the "intensity" value of each pixel is equal 
        to the number of z slices where the binary cell mask is True
    mask_name: str,
        the name or nickname of your mask; this determines how the mask is referred to in the metrics tables
    centering_proj: np.ndarray
        a sum projection of the object you want to use as the center of the distribution where the "intensity" value of each pixel is 
        equal to the number of z slices where the binary nucleus mask is True
    obj_proj: np.ndarray,
        a sum projection of the stuff you want to measure where the "intensity" value of each pixel is equal to the number of z slices 
        where the binary organelle mask is True (for a segmented image) or the total intensity at that point (for a gray scale image)
    obj_name: str,
        the name or nickname of your object being measured; used for labeling columns in the dataframe
    bin_count: int,
        the number of concentric rings, or "bins", to create within the mask
    center_on: bool = False,
        True = distribute the bins from the center of the centering object
        False = distribute the bins from the edge of the centering object
    keep_center_as_bin: bool = True
        True = include the centering object area when creating the bins
        False = do not include the centering object area when creating the bins
    scale: Union[tuple, None]=None
        a tuple of floats representing the real-world dimensions for each image dimension (ZYX)
        

    Measurements
    ------------
    If scale is used, "vox_cnt" is replaced by "vol" and "n_pix_ is replaced by "area" in the titles below.
    If no centering object is provided, the related measurements are omitted
    If no mask object is provided, "mask" is replaced by "img" in the titles below

    object: the nickname of what is being measured (e.g., golgi, golgiXER, ER_img)
    XY_n_bins: number of bins
    XY_bins: list of bin number
    XY_mask_vox_cnt_perbin: number of voxels in the 3D cell mask per bin
    XY_obj_vox_cnt_perbin: number of voxels of the 3D object per bin
    XY_center_vox_cnt_perbin: number of voxels of the 3D centering object per bin
    XY_n_pix_perbin: number of pixels per bin in the XY mask
    XY_portion_pix_perbin: the portion of pixels in the XY mask per bin
    XY_n_wedges: number of wedges
    XY_wedges: list of wedge numbers
    XY_mask_vox_cnt_perwedge: number of voxels in the 3D cell mask per wedge
    XY_obj_vox_cnt_perwedge: number of voxels of the 3D object per wedge
    XY_center_vox_cnt_perwedge: number of voxels of the 3D centering object per wedge
    XY_n_pix_perwedge: number of pixels per wedge in the XY mask
    XY_portion_pix_perwedge: the portion of pixels in the XY mask per bin
    XY_wedges_perbin: list of wedges that have >0 pixels in the mask for all bins
    XY_mask_vox_cnt_wedges_perbin: number of voxels in the 3D cell mask per wedge per bin
    XY_obj_vox_cnt_wedges_perbin:number of voxels of the 3D object per wedge per bin
    XY_center_vox_cnt_wedges_perbin: number of voxels of the 3D centering object per wedge per bin
    XY_n_pix_wedges_perbin: number of pixels per wedge per bin in the XY mask
    XY_mask_cv_perbin: the coefficient of variance of the wedges within each bin for the mask
    XY_obj_cv_perbin: the coefficient of variance of the wedges within each bin for the object segmentation
    XY_center_cv_perbin: the coefficient of variance of the wedges within each bin for the centering object

    
    Returns
    -------------
    tab: (pd.DataFrame) table of measurements of the object distribution
    bin_array: (np.ndarray) mask of the concentric rings to measure distribution from
    wedge_array: (np.ndarray) mask of the wedges (pie slices) that divide each bin into 8 parts
    """
    # other parameters that will stay constant
    nobjects = 1

    # create binary arrays
    center_objects = centering_proj > 0 if centering_proj is not None else None
    mask = (mask_proj>0).astype(np.uint16)


    ################   ################
    ## compute distances and make bins and wedges masks
    ################   ################
    # created normalized distances
    normalized_distance, good_mask, i_center, j_center = get_normalized_distance_and_mask(labels=mask, center_objects=center_objects, center_on=center_on)
    if normalized_distance is None:
        print('normalized_distance returned wrong')

    # create bin mask array
    
    if center_objects is None:
        keep_center_as_bin = True
        center_on = True

    if keep_center_as_bin:
        if center_on:
            bin_array = (normalized_distance * bin_count).astype(int)
        else:
            bin_array= ((normalized_distance * (bin_count-1))+1).astype(int)
            bin_array[center_objects]=0
            bin_array[~good_mask]=0
    else:
        good_mask[center_objects]=0
        if center_on:
            normalized_distance[good_mask] = (normalized_distance[good_mask] - normalized_distance[good_mask].min())/(normalized_distance[good_mask].max() - normalized_distance[good_mask].min())
        bin_array = (normalized_distance * bin_count).astype(int)
            
    bin_array[bin_array > bin_count] = bin_count
    
    # create wedges mask array
    i, j = np.mgrid[0 : mask.shape[0], 0 : mask.shape[1]]
    imask = i[good_mask] > i_center[good_mask]
    jmask = j[good_mask] > j_center[good_mask]
    absmask = abs(i[good_mask] - i_center[good_mask]) > abs(j[good_mask] - j_center[good_mask])
    radial_index = (imask.astype(int) + jmask.astype(int) * 2 + absmask.astype(int) * 4)

    wedge_array = np.zeros_like(good_mask, dtype=int)
    wedge_array[good_mask] = radial_index
    

    ################   ################
    ## get histograms
    ################   ################

    ## These measurements are only using the bins created from the edge of the centering object and including the centering object area
    # number of pixels in the good mask
    ngood_pixels = np.sum(good_mask)

    good_labels = mask[good_mask]

    # whole cell bin and wedge measurements
    mask_arrays = [bin_array, wedge_array]
    sections = [bin_count, 8]
    types = ['bin', 'wedge']

    met_dict = {}

    for array, num, name in zip(mask_arrays, sections, types):
        labels_and_bins = (good_labels - 1, array[good_mask])

        # get count of voxels in each bin from the following images
        met_dict[f"XY_{mask_name}_vox_cnt_per{name}"] = [coo_matrix((mask_proj[good_mask], labels_and_bins), shape=(nobjects, num)).toarray().squeeze().tolist()]
        met_dict[f"XY_obj_vox_cnt_per{name}"] = [coo_matrix((obj_proj[good_mask], labels_and_bins), shape=(nobjects, num)).toarray().squeeze().tolist()]
        # does not create key if condition is false
        if center_objects is not None:
            met_dict[f"XY_center_vox_cnt_per{name}"] = [coo_matrix((centering_proj[good_mask], labels_and_bins), shape=(nobjects, num)).toarray().squeeze().tolist()]

        # same concept, but with an empty array to calculate the number of pixels per bin
        n_pixels = [coo_matrix((np.ones(ngood_pixels), labels_and_bins), (nobjects, num)).toarray().squeeze().tolist()]
        met_dict[f"XY_n_pix_per{name}"] = n_pixels

        # total pixels in the mask
        total_pixels = np.sum(n_pixels, 1)
        total_repeated = np.dstack([total_pixels] * (num))[0]

        # get the proportion of pixels in each bin (*100 to get percentage of cell pixels per bin)
        met_dict[f"XY_portion_pix_per{name}"] = [(n_pixels / total_repeated).squeeze().tolist()]


    # per wedge per bin measurements
    bin_names =[]
    cv_mask = []
    cv_obj = []
    if center_objects is not None:
        cv_center = []
    mask_wedge_perbin = []
    obj_wedge_perbin = []
    if center_objects is not None:
        center_wedge_perbin = []
    pxl_cnt_wedge_perbin = []
    wedges_perbin = []

    for bin in range(bin_count):
        bin_mask = good_mask & (bin_array == bin) # selecting the bin as a mask
        bin_pixels = np.sum(bin_mask) # number of pixels in this bin for downstream calculations

        bin_labels = mask[bin_mask] # selecting portion of the cellmask within this bin

        bin_radial_index = radial_index[bin_array[good_mask] == bin] # selecting the portion of the wedges associated to this bin
        labels_and_radii = (bin_labels - 1, bin_radial_index) # (i,j) for coo_matrix function taking into account the 8 wedges within this bin

        # repeating the calculations above using the wedges instead of the bins
        radial_counts_mask = coo_matrix((mask_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray() # amount of cell mask voxels per wedge in this bin
        radial_counts_obj = coo_matrix((obj_proj[bin_mask], labels_and_radii), (nobjects, 8)).toarray() # amount of object voxels per wedges in this bin
        if center_objects is not None:
            radial_counts_center = coo_matrix((centering_proj[bin_mask], labels_and_radii), (nobjects, 8)).toarray() # amount of centering object voxels per wedges in this bin
        pixel_count = coo_matrix((np.ones(bin_pixels), labels_and_radii), (nobjects, 8)).toarray()

        # safe gaurd against one of the wedges having an area of 0
        # np.ma.masked_array - "Masked values of True exclude the corresponding element from any computation."
        n_mask = pixel_count == 0
        radial_counts = [radial_counts_mask, radial_counts_obj]
        radial_counts += [radial_counts_center] if center_objects is not None else []
        radial_cvs = []
        for count in radial_counts:
            radial_norm = np.ma.masked_array(count / pixel_count, n_mask)
            radial_cv = np.std(radial_norm, 1) / np.mean(radial_norm, 1)
            radial_cv[np.sum(~n_mask, 1) == 0] = 0
            radial_cv.mask = np.sum(~n_mask, 1) == 0
            radial_cvs.append(radial_cv)

        bin_name = bin + 1 if bin > 0 else 1
        wedges_perbin_name = np.ma.masked_array([it+1 for it in range(8)])

        bin_names.append(bin_name)
        cv_mask.append(float(np.mean(radial_cvs[0]))) #convert to float to make importing from csv more straightforward
        cv_obj.append(float(np.mean(radial_cvs[1])))
        if center_objects is not None:
            cv_center.append(float(np.mean(radial_cvs[2])))
        mask_wedge_perbin.append(radial_counts[0].squeeze().tolist())
        obj_wedge_perbin.append(radial_counts[1].squeeze().tolist())
        if center_objects is not None:
            center_wedge_perbin.append(radial_counts[2].squeeze().tolist())
        pxl_cnt_wedge_perbin.append(pixel_count.squeeze().tolist())
        wedges_perbin.append(wedges_perbin_name.data.squeeze().tolist())
    

    ################   ################
    ## create data table and account for scale
    ################   ################
    met_dict_1 = {'object': obj_name,
                    'XY_n_bins': bin_count,
                    'XY_bins': [bin_names]}
    met_dict_2 = dict(list(met_dict.items())[:5])
    met_dict_3 = {'XY_n_wedges': 8,
                    'XY_wedges': str([it+1 for it in range(8)])}
    met_dict_4 = dict(list(met_dict.items())[5:])
    met_dict_5 = {'XY_wedges_perbin': [wedges_perbin],
                f'XY_{mask_name}_vox_cnt_wedges_perbin':[mask_wedge_perbin],
                'XY_obj_vox_cnt_wedges_perbin':[obj_wedge_perbin],
                **({'XY_center_vox_cnt_wedges_perbin': [center_wedge_perbin]} if center_objects is not None else {}),
                'XY_n_pix_wedges_perbin': [pxl_cnt_wedge_perbin],
                f'XY_{mask_name}_cv_perbin':[cv_mask],
                'XY_obj_cv_perbin':[cv_obj],
                **({'XY_center_cv_perbin': [cv_center]} if center_objects is not None else {})}

    dict_combined = dict(itertools.chain(met_dict_1.items(), met_dict_2.items(), met_dict_3.items(), met_dict_4.items(), met_dict_5.items()))
    tab = pd.DataFrame(dict_combined)

    # account for scale
    if scale is not None:
        round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
        tab.insert(loc=0, column="scale", value=f"{round_scale}")
        
        # measurements affected by scale
        vol_mets = [f'XY_{mask_name}_vox_cnt_perbin', 'XY_obj_vox_cnt_perbin',
                    *(['XY_center_vox_cnt_perbin'] if center_objects is not None else []), f'XY_{mask_name}_vox_cnt_perwedge',
                    'XY_obj_vox_cnt_perwedge', *(['XY_center_vox_cnt_perwedge'] if center_objects is not None else []),
                    f'XY_{mask_name}_vox_cnt_wedges_perbin', 'XY_obj_vox_cnt_wedges_perbin',
                    *(['XY_center_vox_cnt_wedges_perbin'] if center_objects is not None else [])]
        
        area_mets = ['XY_n_pix_perbin', 'XY_n_pix_perwedge', 'XY_n_pix_wedges_perbin']

        for met in vol_mets:
            tab[met.replace('_vox_cnt_', "_vol_")] = [(np.float_(tab[met][0]) * np.prod(scale)).squeeze().tolist()]
        for met in area_mets:
            tab[met.replace('_n_pix_', "_area_")] = [(np.float_(tab[met][0]) * np.prod(scale[1:])).squeeze().tolist()]

    else: 
        tab.insert(loc=0, column="scale", value=f"{tuple(np.ones(3))}")
    # add mask name to table
    tab.insert(loc=0, column = "mask_name", value = mask_name)
    
    return tab, bin_array, wedge_array

# def get_concentric_distribution(
#         mask_proj: np.ndarray,
#         centering_proj: np.ndarray,
#         obj_proj: np.ndarray,
#         obj_name: str,
#         bin_count: int,
#         center_on: bool = False,
#         keep_center_as_bin: bool = True,
#         scale: Union[tuple, None]=None):
#     """
#     Based on CellProfiler's measureobjectintensitydistribution. Measure the distribution of segmented objects within a masked area. 
#     In our case, we will usually utilize this function to measure the amount of an organelle within the cell.
#     Radial bins are created out from a center point, usually the nucleus edge.

    
#     Parameters
#     ------------
#     mask_proj: np.ndarray
#         a sum projection of the region you want to measure the distribution from where the "intensity" value of each pixel is equal 
#         to the number of z slices where the binary cell mask is True
#     centering_proj: np.ndarray
#         a sum projection of the object you want to use as the center of the distribution where the "intensity" value of each pixel is 
#         equal to the number of z slices where the binary nucleus mask is True
#     obj_proj: np.ndarray,
#         a sum projection of the stuff you want to measure where the "intensity" value of each pixel is equal to the number of z slices 
#         where the binary organelle mask is True (for a segmented image) or the total intensity at that point (for a gray scale image)
#     obj_name: str,
#         the name or nickname of your object being measured; used for labeling columns in the dataframe
#     bin_count: int,
#         the number of concentric rings, or "bins", to create within the mask
#     center_on: bool = False,
#         True = distribute the bins from the center of the centering object
#         False = distribute the bins from the edge of the centering object
#     keep_center_as_bin: bool = True
#         True = include the centering object area when creating the bins
#         False = do not include the centering object area when creating the bins
#     scale: Union[tuple, None]=None
#         a tuple of floats representing the real-world dimensions for each image dimension (ZYX)
        

#     Measurements
#     ------------
#     If scale is used, "vox_cnt" is replaced by "vol" and "n_pix_ is replaced by "area" in the title below.

#     object: the nickname of what is being measured (e.g., golgi, golgiXER, ER_img)
#     XY_n_bins: number of bins
#     XY_bins: list of bin number
#     XY_mask_vox_cnt_perbin: number of voxels in the 3D cell mask per bin
#     XY_obj_vox_cnt_perbin: number of voxels of the 3D object per bin
#     XY_center_vox_cnt_perbin: number of voxels of the 3D centering object per bin
#     XY_n_pix_perbin: number of pixels per bin in the XY mask
#     XY_portion_pix_perbin: the portion of pixels in the XY mask per bin
#     XY_n_wedges: number of wedges
#     XY_wedges: list of wedge numbers
#     XY_mask_vox_cnt_perwedge: number of voxels in the 3D cell mask per wedge
#     XY_obj_vox_cnt_perwedge: number of voxels of the 3D object per wedge
#     XY_center_vox_cnt_perwedge: number of voxels of the 3D centering object per wedge
#     XY_n_pix_perwedge: number of pixels per wedge in the XY mask
#     XY_portion_pix_perwedge: the portion of pixels in the XY mask per bin
#     XY_wedges_perbin: list of wedges that have >0 pixels in the mask for all bins
#     XY_mask_vox_cnt_wedges_perbin: number of voxels in the 3D cell mask per wedge per bin
#     XY_obj_vox_cnt_wedges_perbin:number of voxels of the 3D object per wedge per bin
#     XY_center_vox_cnt_wedges_perbin: number of voxels of the 3D centering object per wedge per bin
#     XY_n_pix_wedges_perbin: number of pixels per wedge per bin in the XY mask
#     XY_mask_cv_perbin: the coefficient of variance of the wedges within each bin for the mask
#     XY_obj_cv_perbin: the coefficient of variance of the wedges within each bin for the object segmentation
#     XY_center_cv_perbin: the coefficient of variance of the wedges within each bin for the centering object

    
#     Returns
#     -------------
#     tab: (pd.DataFrame) table of measurements of the object distribution
#     bin_array: (np.ndarray) mask of the concentric rings to measure distribution from
#     wedge_array: (np.ndarray) mask of the wedges (pie slices) that divide each bin into 8 parts
#     """
#     # other parameters that will stay constant
#     nobjects = 1

#     # create binary arrays
#     center_objects = centering_proj > 0 if centering_proj is not None else None
#     mask = (mask_proj>0).astype(np.uint16)


#     ################   ################
#     ## compute distances and make bins and wedges masks
#     ################   ################
#     # created normalized distances
#     normalized_distance, good_mask, i_center, j_center = get_normalized_distance_and_mask(labels=mask, center_objects=center_objects, center_on=center_on)
#     if normalized_distance is None:
#         print('normalized_distance returned wrong')

#     # create bin mask array
    
#     if center_objects is None:
#         keep_center_as_bin = True
#         center_on = True

#     if keep_center_as_bin:
#         if center_on:
#             bin_array = (normalized_distance * bin_count).astype(int)
#         else:
#             bin_array= ((normalized_distance * (bin_count-1))+1).astype(int)
#             bin_array[center_objects]=0
#             bin_array[~good_mask]=0
#     else:
#         good_mask[center_objects]=0
#         if center_on:
#             normalized_distance[good_mask] = (normalized_distance[good_mask] - normalized_distance[good_mask].min())/(normalized_distance[good_mask].max() - normalized_distance[good_mask].min())
#         bin_array = (normalized_distance * bin_count).astype(int)
            
#     bin_array[bin_array > bin_count] = bin_count
    
#     # create wedges mask array
#     i, j = np.mgrid[0 : mask.shape[0], 0 : mask.shape[1]]
#     imask = i[good_mask] > i_center[good_mask]
#     jmask = j[good_mask] > j_center[good_mask]
#     absmask = abs(i[good_mask] - i_center[good_mask]) > abs(j[good_mask] - j_center[good_mask])
#     radial_index = (imask.astype(int) + jmask.astype(int) * 2 + absmask.astype(int) * 4)

#     wedge_array = np.zeros_like(good_mask, dtype=int)
#     wedge_array[good_mask] = radial_index
    

#     ################   ################
#     ## get histograms
#     ################   ################

#     # number of pixels in the good mask
#     ngood_pixels = np.sum(good_mask)

#     good_labels = mask[good_mask]

#     # whole cell bin and wedge measurements
#     mask_arrays = [bin_array, wedge_array]
#     sections = [bin_count, 8]
#     types = ['bin', 'wedge']

#     met_dict = {}

#     for array, num, name in zip(mask_arrays, sections, types):
#         labels_and_bins = (good_labels - 1, array[good_mask])

#         # get count of voxels in each bin from the following images
#         met_dict[f"XY_mask_vox_cnt_per{name}"] = [coo_matrix((mask_proj[good_mask], labels_and_bins), shape=(nobjects, num)).toarray().squeeze().tolist()]
#         met_dict[f"XY_obj_vox_cnt_per{name}"] = [coo_matrix((obj_proj[good_mask], labels_and_bins), shape=(nobjects, num)).toarray().squeeze().tolist()]
#         # does not create key if condition is false
#         if center_objects is not None:
#             met_dict[f"XY_center_vox_cnt_per{name}"] = [coo_matrix((centering_proj[good_mask], labels_and_bins), shape=(nobjects, num)).toarray().squeeze().tolist()]

#         # same concept, but with an empty array to calculate the number of pixels per bin
#         n_pixels = [coo_matrix((np.ones(ngood_pixels), labels_and_bins), (nobjects, num)).toarray().squeeze().tolist()]
#         met_dict[f"XY_n_pix_per{name}"] = n_pixels

#         # total pixels in the mask
#         total_pixels = np.sum(n_pixels, 1)
#         total_repeated = np.dstack([total_pixels] * (num))[0]

#         # get the proportion of pixels in each bin (*100 to get percentage of cell pixels per bin)
#         met_dict[f"XY_portion_pix_per{name}"] = [(n_pixels / total_repeated).squeeze().tolist()]


#     # per wedge per bin measurements
#     bin_names =[]
#     cv_mask = []
#     cv_obj = []
#     if center_objects is not None:
#         cv_center = []
#     mask_wedge_perbin = []
#     obj_wedge_perbin = []
#     if center_objects is not None:
#         center_wedge_perbin = []
#     pxl_cnt_wedge_perbin = []
#     wedges_perbin = []

#     for bin in range(bin_count):
#         bin_mask = good_mask & (bin_array == bin) # selecting the bin as a mask
#         bin_pixels = np.sum(bin_mask) # number of pixels in this bin for downstream calculations

#         bin_labels = mask[bin_mask] # selecting portion of the cellmask within this bin

#         bin_radial_index = radial_index[bin_array[good_mask] == bin] # selecting the portion of the wedges associated to this bin
#         labels_and_radii = (bin_labels - 1, bin_radial_index) # (i,j) for coo_matrix function taking into account the 8 wedges within this bin

#         # repeating the calculations above using the wedges instead of the bins
#         radial_counts_mask = coo_matrix((mask_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray() # amount of cell mask voxels per wedge in this bin
#         radial_counts_obj = coo_matrix((obj_proj[bin_mask], labels_and_radii), (nobjects, 8)).toarray() # amount of object voxels per wedges in this bin
#         if center_objects is not None:
#             radial_counts_center = coo_matrix((centering_proj[bin_mask], labels_and_radii), (nobjects, 8)).toarray() # amount of centering object voxels per wedges in this bin
#         pixel_count = coo_matrix((np.ones(bin_pixels), labels_and_radii), (nobjects, 8)).toarray()

#         # safe gaurd against one of the wedges having an area of 0
#         # np.ma.masked_array - "Masked values of True exclude the corresponding element from any computation."
#         n_mask = pixel_count == 0
#         radial_counts = [radial_counts_mask, radial_counts_obj]
#         radial_counts += [radial_counts_center] if center_objects is not None else []
#         radial_cvs = []
#         for count in radial_counts:
#             radial_norm = np.ma.masked_array(count / pixel_count, n_mask)
#             radial_cv = np.std(radial_norm, 1) / np.mean(radial_norm, 1)
#             radial_cv[np.sum(~n_mask, 1) == 0] = 0
#             radial_cv.mask = np.sum(~n_mask, 1) == 0
#             radial_cvs.append(radial_cv)

#         bin_name = bin + 1 if bin > 0 else 1
#         wedges_perbin_name = np.ma.masked_array([it+1 for it in range(8)])

#         bin_names.append(bin_name)
#         cv_mask.append(float(np.mean(radial_cvs[0]))) #convert to float to make importing from csv more straightforward
#         cv_obj.append(float(np.mean(radial_cvs[1])))
#         if center_objects is not None:
#             cv_center.append(float(np.mean(radial_cvs[2])))
#         mask_wedge_perbin.append(radial_counts[0].squeeze().tolist())
#         obj_wedge_perbin.append(radial_counts[1].squeeze().tolist())
#         if center_objects is not None:
#             center_wedge_perbin.append(radial_counts[2].squeeze().tolist())
#         pxl_cnt_wedge_perbin.append(pixel_count.squeeze().tolist())
#         wedges_perbin.append(wedges_perbin_name.data.squeeze().tolist())
    

#     ################   ################
#     ## create data table and account for scale
#     ################   ################
#     met_dict_1 = {'object': obj_name,
#                     'XY_n_bins': bin_count,
#                     'XY_bins': [bin_names]}
#     met_dict_2 = dict(list(met_dict.items())[:5])
#     met_dict_3 = {'XY_n_wedges': 8,
#                     'XY_wedges': str([it+1 for it in range(8)])}
#     met_dict_4 = dict(list(met_dict.items())[5:])
#     met_dict_5 = {'XY_wedges_perbin': [wedges_perbin],
#                 'XY_mask_vox_cnt_wedges_perbin':[mask_wedge_perbin],
#                 'XY_obj_vox_cnt_wedges_perbin':[obj_wedge_perbin],
#                 **({'XY_center_vox_cnt_wedges_perbin': [center_wedge_perbin]} if center_objects is not None else {}),
#                 'XY_n_pix_wedges_perbin': [pxl_cnt_wedge_perbin],
#                 'XY_mask_cv_perbin':[cv_mask],
#                 'XY_obj_cv_perbin':[cv_obj],
#                 **({'XY_center_cv_perbin': [cv_center]} if center_objects is not None else {})}

#     dict_combined = dict(itertools.chain(met_dict_1.items(), met_dict_2.items(), met_dict_3.items(), met_dict_4.items(), met_dict_5.items()))
#     tab = pd.DataFrame(dict_combined)

#     # account for scale
#     if scale is not None:
#         round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
#         tab.insert(loc=1, column="scale", value=f"{round_scale}")
        
#         # measurements affected by scale
#         vol_mets = ['XY_mask_vox_cnt_perbin', 'XY_obj_vox_cnt_perbin',
#                     *(['XY_center_vox_cnt_perbin'] if center_objects is not None else []), 'XY_mask_vox_cnt_perwedge',
#                     'XY_obj_vox_cnt_perwedge', *(['XY_center_vox_cnt_perwedge'] if center_objects is not None else []),
#                     'XY_mask_vox_cnt_wedges_perbin', 'XY_obj_vox_cnt_wedges_perbin',
#                     *(['XY_center_vox_cnt_wedges_perbin'] if center_objects is not None else [])]
        
#         area_mets = ['XY_n_pix_perbin', 'XY_n_pix_perwedge', 'XY_n_pix_wedges_perbin']

#         for met in vol_mets:
#             tab[met.replace('_vox_cnt_', "_vol_")] = [(np.float_(tab[met][0]) * np.prod(scale)).squeeze().tolist()]
#         for met in area_mets:
#             tab[met.replace('_n_pix_', "_area_")] = [(np.float_(tab[met][0]) * np.prod(scale[1:])).squeeze().tolist()]

#     else: 
#         tab.insert(loc=2, column="scale", value=f"{tuple(np.ones(3))}")

#     return tab, bin_array, wedge_array

# Zernicke routines.  inspired by cellprofiler, but heavily simplified
### USED ###
def zernike_metrics(pixels,z):
    """
    
    """
    vr = np.sum(pixels[:,:,np.newaxis]*z.real, axis=(0,1))
    vi = np.sum(pixels[:,:,np.newaxis]*z.imag, axis=(0,1))    
    magnitude = np.sqrt(vr * vr + vi * vi) / pixels.sum()
    phase = np.arctan2(vr, vi)
    # return {"zer_mag": magnitude, "zer_phs": phase}
    return magnitude, phase


## USED ###
def zernike_polynomial(labels, zernike_is):
    """
    

    """
    # First, get a table of centers and radii of minimum enclosing
    # circles for the cellmask
    ij, r = centrosome.cpmorphology.minimum_enclosing_circle( labels )
    # Then compute x and y, the position of each labeled pixel
    # within a unit circle around the object
    iii, jjj = np.mgrid[0 : labels.shape[0], 0 : labels.shape[1]]

    # translate+scale
    iii = (iii-ij[0][0] ) / r
    jjj = (jjj-ij[0][1] ) / r

    z = centrosome.zernike.construct_zernike_polynomials(
        iii, jjj, zernike_is
    )
    return z
    


### USED ###
def get_zernike_metrics(        
        mask_proj: np.ndarray,
        mask_name: str,
        centering_proj: Union[np.ndarray, None], 
        obj_proj: np.ndarray,
        obj_name: str,
        zernike_degree: int = 9 ):

    """
    Compute Zernike-based metrics for a projected mask and corresponding object projection.
    This function computes Zernike polynomials over regions derived from the
    ``mask_proj`` and then measures the Zernike magnitudes and phases for:
    - the mask projection (``mask_proj``),
    - the object projection (``obj_proj``),
    - the centering projection (``centering_proj``), if provided.
    The results are returned as a single-row :class:`pandas.DataFrame` that contains
    the Zernike indices (``n`` and ``m``) together with the magnitude and phase
    vectors for each input.
    
    Parameters
    ----------
    mask_proj : numpy.ndarray
        a sum projection of the region you want to measure the distribution from where the "intensity" value of each pixel is equal 
        to the number of z slices where the binary cell mask is True
    mask_name : str
        the name or nickname of your mask; this determines how the mask is referred to in the metrics tables
    centering_proj : numpy.ndarray or None
        Optional 2D array used for centering. If not ``None``, Zernike metrics
        are also computed for this projection and added to the output.
    obj_proj : numpy.ndarray
        2D array representing the projection of the object of interest
        for which Zernike metrics will be computed using the same Zernike basis.
    obj_name : str
        Name of the object or channel represented by ``obj_proj``; stored in the
        ``"object"`` column of the output DataFrame.
    zernike_degree : int, optional
        Maximum degree of the Zernike polynomials. All Zernike indices up to
        and including this degree (plus one in the underlying library call) are
        used. Defaults to 9.
    Returns
    -------
    pandas.DataFrame
        A single-row DataFrame with Zernike information. Columns include:
        - ``"object"``: the provided ``obj_name``.
        - ``"zernike_n"`` and ``"zernike_m"``: lists of the Zernike index pairs.
        - ``f"zernike_{mask_name}_mag"`` and ``f"zernike_{mask_name}_phs"``:
          lists of magnitudes and phases computed from ``mask_proj``.
        - ``"zernike_obj_mag"`` and ``"zernike_obj_phs"``: magnitudes and phases
          for ``obj_proj``.
        - ``"zernike_center_mag"`` and ``"zernike_center_phs"``: magnitudes and
          phases for ``centering_proj``, if ``centering_proj`` is not ``None``.
    
    """
    

    labels = label(mask_proj>0) #extent as 0,1 rather than bool
    zernike_indexes = centrosome.zernike.get_zernike_indexes( zernike_degree + 1)


    z = zernike_polynomial(labels, zernike_indexes)

    z_m = zernike_metrics(mask_proj, z)
    z_obj = zernike_metrics(obj_proj, z)
    if centering_proj is not None:
        z_c = zernike_metrics(centering_proj, z)


    # nm_labels = [f"{n}_{m}" for (n, m) in (zernike_indexes)
    stats_tab = pd.DataFrame({'object':obj_name,
                                'zernike_n':[zernike_indexes[:,0].tolist()],
                                'zernike_m':[zernike_indexes[:,1].tolist()],
                                f'zernike_{mask_name}_mag':[z_m[0].tolist()],
                                f'zernike_{mask_name}_phs':[z_m[1].tolist()],   
                                'zernike_obj_mag':[z_obj[0].tolist()],
                                'zernike_obj_phs':[z_obj[1].tolist()],
                                **({'zernike_center_mag':[z_c[0].tolist()]} if centering_proj is not None else {}),
                                **({'zernike_center_phs':[z_c[1].tolist()]} if centering_proj is not None else {})})

    return stats_tab

# def get_zernike_metrics(        
#         cellmask_proj: np.ndarray,
#         nucleus_proj: Union[np.ndarray, None], 
#         org_proj: np.ndarray,
#         organelle_name: str,
#         zernike_degree: int = 9 ):

#     """
    
#     """

#     labels = label(cellmask_proj>0) #extent as 0,1 rather than bool
#     zernike_indexes = centrosome.zernike.get_zernike_indexes( zernike_degree + 1)


#     z = zernike_polynomial(labels, zernike_indexes)

#     z_cm = zernike_metrics(cellmask_proj, z)
#     z_org = zernike_metrics(org_proj, z)
#     if nucleus_proj is not None:
#         z_nuc = zernike_metrics(nucleus_proj, z)



#     # nm_labels = [f"{n}_{m}" for (n, m) in (zernike_indexes)
#     stats_tab = pd.DataFrame({'object':organelle_name,
#                                 'zernike_n':[zernike_indexes[:,0].tolist()],
#                                 'zernike_m':[zernike_indexes[:,1].tolist()],
#                                 'zernike_mask_mag':[z_cm[0].tolist()],
#                                 'zernike_mask_phs':[z_cm[1].tolist()],   
#                                 'zernike_obj_mag':[z_org[0].tolist()],
#                                 'zernike_obj_phs':[z_org[1].tolist()],
#                                 **({'zernike_center_mag':[z_nuc[0].tolist()]} if nucleus_proj is not None else {}),
#                                 **({'zernike_center_phs':[z_nuc[1].tolist()]} if nucleus_proj is not None else {})})

#     return stats_tab



### USED ###
def get_XY_distribution(        
        mask: Union[np.ndarray,None],
        mask_name: str,
        centering_obj: Union[np.ndarray,None],
        obj:np.ndarray,
        obj_name: str,
        scale: Union[tuple, None]=None,
        num_bins: Union[int, None] = 5,
        center_on: bool = False,
        keep_center_as_bin: bool = True,
        zernike_degrees: Union[int, None] = None):

    """
    Params
    ----------
    mask: np.ndarray,
        a binary 3D (ZYX) np.ndarray of the area that will be measured from
    mask_name: str
        the name or nickname for the mask object; this name will appear in the metrics output
    centering_obj: np.ndarray
        a binary 3D (ZYX) np.ndarray of the object that will be used as the center of the concentric rins ("bins")
    obj: np.ndarray
        a 3D (ZYX) np.ndarray image of what will be measured within the masked area
    obj_name: str
        the name or nickname for the obj being measured; this will appear as a column in the output datasheet
    scale: Union[tuple, None]=None
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)
    num_bins: Union[int,None] = None
        the number of concentric rings to draw between the centering object and edge of the mask; None will result in 5 bins
    center_on: bool = False
        True = distribute the bins from the center of the centering object
        False = distribute the bins from the edge of the centering object
    keep_center_as_bin: bool = True
        True = include the centering object area when creating the bins
        False = do not include the centering object area when creating the bins
    zernike_degrees: Union[int,None] = None
        the number of Zernike degrees to include for the Zernike shape descriptors; if None, the Zernike measurements will not 
        be included in the output


    Returns
    -----------
    XY_metrics:
        a pandas Dataframe of bin, wedge, and Zernike measurements
    dist_bin_mask:
        an np.ndarray mask of the concentric ring bins
    dist_wedge_mask 
        an np.ndarray mask of the 8 radial wedges

    """
    # create sum Z projections
    # the mask that will be applied to the centering and organelle object
    m = mask.astype(bool) if mask is not None else None

    center_proj = create_masked_sum_projection(centering_obj,m) if centering_obj is not None else None
    obj_proj = create_masked_sum_projection(obj,m)

    # mask 2d sum projection
    mask_proj = create_masked_sum_projection(mask) if mask is not None else np.full_like(obj_proj,obj.shape[0])
 

    XY_metrics, dist_bin_mask, dist_wedge_mask = get_concentric_distribution(mask_proj=mask_proj,
                                                        mask_name = mask_name,                     
                                                        centering_proj=center_proj, 
                                                        obj_proj=obj_proj, 
                                                        obj_name=obj_name, 
                                                        scale=scale,
                                                        bin_count=num_bins, 
                                                        center_on=center_on,
                                                        keep_center_as_bin=keep_center_as_bin)
    
    if zernike_degrees is not None:
        zernike_metrics = get_zernike_metrics(mask_proj=mask_proj,
                                            mask_name = mask_name, 
                                            obj_proj=obj_proj,
                                            obj_name=obj_name, 
                                            centering_proj=center_proj, 
                                            zernike_degree=zernike_degrees)
        
        XY_metrics = pd.merge(XY_metrics, zernike_metrics, on="object")

    return XY_metrics, dist_bin_mask, dist_wedge_mask

# def get_XY_distribution(        
#         mask: np.ndarray,
#         centering_obj: np.ndarray,
#         obj:np.ndarray,
#         obj_name: str,
#         scale: Union[tuple, None]=None,
#         num_bins: Union[int, None] = 5,
#         center_on: bool = False,
#         keep_center_as_bin: bool = True,
#         zernike_degrees: Union[int, None] = None):

#     """
#     Params
#     ----------
#     mask_obj: np.ndarray,
#         a binary 3D (ZYX) np.ndarray of the area that will be measured from
#     centering_obj: np.ndarray
#         a binary 3D (ZYX) np.ndarray of the object that will be used as the center of the concentric rins ("bins")
#     obj: np.ndarray
#         a 3D (ZYX) np.ndarray image of what will be measured within the masked area
#     obj_name: str
#         the name or nickname for the obj being measured; this will appear as a column in the output datasheet
#     scale: Union[tuple, None]=None
#         a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)
#     num_bins: Union[int,None] = None
#         the number of concentric rings to draw between the centering object and edge of the mask; None will result in 5 bins
#     center_on: bool = False
#         True = distribute the bins from the center of the centering object
#         False = distribute the bins from the edge of the centering object
#     keep_center_as_bin: bool = True
#         True = include the centering object area when creating the bins
#         False = do not include the centering object area when creating the bins
#     zernike_degrees: Union[int,None] = None
#         the number of zernike degrees to include for the zernike shape descriptors; if None, the zernike measurements will not 
#         be included in the output


#     Returns
#     -----------
#     XY_metrics:
#         a pandas Dataframe of bin, wedge, and zernike measurements
#     dist_bin_mask:
#         an np.ndarray mask of the concentric ring bins
#     dist_wedge_mask 
#         an np.ndarray mask of the 8 radial wedges

#     """

#     mask_proj = create_masked_sum_projection(mask)
#     center_proj = create_masked_sum_projection(centering_obj,mask.astype(bool)) if centering_obj is not None else None
#     obj_proj = create_masked_sum_projection(obj,mask.astype(bool))
 

#     XY_metrics, dist_bin_mask, dist_wedge_mask = get_concentric_distribution(mask_proj=mask_proj, 
#                                                         centering_proj=center_proj, 
#                                                         obj_proj=obj_proj, 
#                                                         obj_name=obj_name, 
#                                                         scale=scale,
#                                                         bin_count=num_bins, 
#                                                         center_on=center_on,
#                                                         keep_center_as_bin=keep_center_as_bin)
    
#     if zernike_degrees is not None:
#         zernike_metrics = get_zernike_metrics(cellmask_proj=mask_proj, 
#                                             org_proj=obj_proj,
#                                             organelle_name=obj_name, 
#                                             nucleus_proj=center_proj, 
#                                             zernike_degree=zernike_degrees)
        
#         XY_metrics = pd.merge(XY_metrics, zernike_metrics, on="object")

#     return XY_metrics, dist_bin_mask, dist_wedge_mask

###################################
### Z DISTRIBUTION
###################################

### USED ###
def create_masked_depth_projection(img_in:np.ndarray, mask:Union[np.ndarray, None]=None, to_bool:bool=True) -> np.ndarray:
    """
    create a masked projection by summing together all XY pixels per Z plane/slice
    """
    img_out = img_in.astype(bool) if to_bool else img_in
    if mask is not None:
        img_out = apply_mask(img_out, mask)
    
    return img_out.sum(axis=(1,2))

### USED ###
def get_Z_distribution(        
        mask: Union[np.ndarray,None],
        mask_name: str,
        obj:np.ndarray,
        obj_name: str,
        center_obj: Union[np.ndarray, None],
        scale: Union[tuple, None] = None
        ):
    """
    quantification of distribution along the Z axis; all XY pixels are summed together per Z slice and then quantified

    Parameters
    ------------
    mask_obj: np.ndarray,
        a binary 3D (ZYX) np.ndarray of the area that will be measured from
    mask_name: str
        the name or nickname for the mask object; this name will appear in the metrics output
    obj: np.ndarray
        a 3D (ZYX) np.ndarray image of what will be measured within the masked area
    obj_name: str
        the name or nickname for the obj being measured; this will appear as a column in the output datasheet
    centering_obj: np.ndarray
        optional - a binary 3D (ZYX) np.ndarray utilized as the center/reference point of the area; for cells, this is usually the nucleus
    scale: Union[tuple, None]=None
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)

    Returns
    -----------
    Z_tab:
        a pandas Dataframe of measurements for each z slice

    """
    # the mask that will be applied to the centering and organelle object
    m = mask.astype(bool) if mask is not None else None

    # flattened
    obj_proj = create_masked_depth_projection(obj, m)
    mask_proj = create_masked_depth_projection(mask) if mask is not None else np.full_like(obj_proj, np.prod(obj.shape[1:]))
    center_proj = create_masked_depth_projection(center_obj, m) if center_obj is not None else None

    Zdist_tab = pd.DataFrame({'object':obj_name,
                            # non-scaled measurments
                            'Z_n_slices':obj.shape[0],
                            'Z_slices':[[i for i in range(obj.shape[0])]],
                            f'Z_{mask_name}_vox_cnt':[mask_proj.tolist()],
                            'Z_obj_vox_cnt':[obj_proj.tolist()],
                            **({'Z_center_vox_cnt': [center_proj.tolist()]} if center_proj is not None else {})
                        })
    
    # scaled measurements added if applicable
    if scale is not None:
        round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
        Zdist_tab.insert(loc=0, column="scale", value=f"{round_scale}")

        Zdist_tab['Z_height'] = obj.shape[0] * scale[0]
        Zdist_tab[f'Z_{mask_name}_volume'] = [(mask_proj * np.prod(scale)).tolist()]
        Zdist_tab['Z_obj_volume'] = [(obj_proj * np.prod(scale)).tolist()]
        if center_proj is not None:
            Zdist_tab['Z_center_volume'] = [(center_proj * np.prod(scale)).tolist()]
    else: 
        Zdist_tab.insert(loc=0, column="scale", value=f"{tuple(np.ones(3))}")
    Zdist_tab.insert(loc=0, column = "mask_name", value = mask_name)
    return Zdist_tab

# def get_Z_distribution(        
#         mask: np.ndarray,
#         obj:np.ndarray,
#         obj_name: str,
#         center_obj: Union[np.ndarray, None],
#         scale: Union[tuple, None] = None
#         ):
#     """
#     quantification of distribution along the Z axis; all XY pixels are summed together per Z slice and then quantified

#     Parameters
#     ------------
#     mask_obj: np.ndarray,
#         a binary 3D (ZYX) np.ndarray of the area that will be measured from
#     obj: np.ndarray
#         a 3D (ZYX) np.ndarray image of what will be measured within the masked area
#     obj_name: str
#         the name or nickname for the obj being measured; this will appear as a column in the output datasheet
#     centering_obj: np.ndarray
#         optional - a binary 3D (ZYX) np.ndarray utilized as the center/reference point of the area; for cells, this is usually the nucleus
#     scale: Union[tuple, None]=None
#         a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)

#     Returns
#     -----------
#     Z_tab:
#         a pandas Dataframe of measurements for each z slice

#     """

#     # flattened
#     mask_proj = create_masked_depth_projection(mask)
#     obj_proj = create_masked_depth_projection(obj, mask.astype(bool))
#     center_proj = create_masked_depth_projection(center_obj, mask.astype(bool)) if center_obj is not None else None

#     Zdist_tab = pd.DataFrame({'object':obj_name,
#                             # non-scaled measurments
#                             'Z_n_slices':mask.shape[0],
#                             'Z_slices':[[i for i in range(mask.shape[0])]],
#                             'Z_mask_vox_cnt':[mask_proj.tolist()],
#                             'Z_obj_vox_cnt':[obj_proj.tolist()]})
#     if center_proj is not None:
#         Zdist_tab['Z_center_vox_cnt'] = [center_proj.tolist()]
    
#     # scaled measurements added if applicable
#     if scale is not None:
#         round_scale = (round(scale[0], 4), round(scale[1], 4), round(scale[2], 4))
#         Zdist_tab.insert(loc=1, column="scale", value=f"{round_scale}")

#         Zdist_tab['Z_height'] = mask.shape[0] * scale[0]
#         Zdist_tab['Z_mask_volume'] = [(mask_proj * np.prod(scale)).tolist()]
#         Zdist_tab['Z_obj_volume'] = [(obj_proj * np.prod(scale)).tolist()]
#         if center_proj is not None:
#             Zdist_tab['Z_center_volume'] = [(center_proj * np.prod(scale)).tolist()]
#     else: 
#         Zdist_tab.insert(loc=2, column="scale", value=f"{tuple(np.ones(3))}")

#     return Zdist_tab


########## DEPRICATED ##########
# def get_normalized_distance_and_mask(labels: np.ndarray, 
#                                       center_objects: Union[np.ndarray, None], 
#                                       center_on: bool):
#     """
#     helper for radial distribution
#     Parameters:
#     ----------
#     labels:
#         2D (YX) np.ndarray - normally the result of a binary ZYX segmentation of the cell mask after a sum projection across the Z dimension
#     center_object:
#         2D (YX) np.ndarray - normally the result of a binary ZYX segmentation of the nucleus after a sum projection across the Z dimension.
#         If no centering object is included, the center of the labels will be used.
#     center_on:
#         True = the center of the centering object will be used as the starting point to calculate the distance from the center
#         False = the edge of the centering object will be used as the starting point to calculate the distance from the center
    
#     Output:
#     ----------
#     normalized_distance:
#         2D (YX) np.ndarray with intensity values representing the distance btween the edge of the "labels" and the centering object
#     good_mask:
#         mask of the areas that were included in the normalized_distance output
#     i_center
#     j_center
#     """

#     d_to_edge = centrosome.cpmorphology.distance_to_edge(labels)

#     if center_objects is not None:
#         center_labels = label(center_objects)
#         pixel_counts = centrosome.cpmorphology.fixup_scipy_ndimage_result(ndi_sum(np.ones(center_labels.shape), 
#                                                                                   center_labels, 
#                                                                                   np.arange(1, np.max(center_labels) + 1, dtype=np.int32)))
#         good = pixel_counts > 0
#         i, j = (centrosome.cpmorphology.centers_of_labels(center_labels) + 0.5).astype(int)
#         ig = i[good]
#         jg = j[good]
#         lg = np.arange(1, len(i) + 1)[good]
        
#         if center_on:  # Reduce the propagation labels to the centers of the centering objects
#             center_labels = np.zeros(center_labels.shape, int)
#             center_labels[ig, jg] = lg

#         cl, d_from_center = centrosome.propagate.propagate(np.zeros(center_labels.shape), center_labels, labels != 0, 1)
#         cl[labels == 0] = 0

#         missing_mask = (labels != 0) & (cl == 0)
#         missing_labels = np.unique(labels[missing_mask])
        
#         if len(missing_labels):
#             print("WTF!!  how did we have missing labels?")
#             all_centers = centrosome.cpmorphology.centers_of_labels(labels)
#             missing_i_centers, missing_j_centers = all_centers[:, missing_labels-1]
#             di = missing_i_centers[:, np.newaxis] - ig[np.newaxis, :]
#             dj = missing_j_centers[:, np.newaxis] - jg[np.newaxis, :]
#             missing_best = lg[np.argsort(di * di + dj * dj)[:, 0]]
#             best = np.zeros(np.max(labels) + 1, int)
#             best[missing_labels] = missing_best
#             cl[missing_mask] = best[labels[missing_mask]]

#             iii, jjj = np.mgrid[0 : labels.shape[0], 0 : labels.shape[1]]
#             di = iii[missing_mask] - i[cl[missing_mask] - 1]
#             dj = jjj[missing_mask] - j[cl[missing_mask] - 1]
#             d_from_center[missing_mask] = np.sqrt(di * di + dj * dj)

#         good_mask = cl > 0
            
#     else:
#         i, j = centrosome.cpmorphology.maximum_position_of_labels(d_to_edge, labels, [1])
#         center_labels = np.zeros(labels.shape, int)
#         center_labels[i, j] = labels[i, j]
#         colors = centrosome.cpmorphology.color_labels(labels)
#         ncolors = np.max(colors)
#         d_from_center = np.zeros(labels.shape)
#         cl = np.zeros(labels.shape, int)

#         for color in range(1, ncolors + 1):
#             mask = colors == color
#             l, d = centrosome.propagate.propagate( np.zeros(center_labels.shape), center_labels, mask, 1)
#             d_from_center[mask] = d[mask]
#             cl[mask] = l[mask]

#         good_mask = cl > 0

#     i_center = np.zeros(cl.shape)
#     i_center[good_mask] = i[cl[good_mask] - 1]

#     j_center = np.zeros(cl.shape)
#     j_center[good_mask] = j[cl[good_mask] - 1]

#     normalized_distance = np.zeros(labels.shape)
#     total_distance = d_from_center + d_to_edge
#     normalized_distance[good_mask] = d_from_center[good_mask] / (total_distance[good_mask] + 0.001)
    
#     return normalized_distance, good_mask, i_center, j_center
   
# def get_zernike_metrics(        
#         cellmask_proj: np.ndarray,
#         nucleus_proj: Union[np.ndarray, None], 
#         org_proj: np.ndarray,
#         organelle_name: str,
#         zernike_degree: int = 9 
#         ):

#     """
    
#     """

#     labels = label(cellmask_proj>0) #extent as 0,1 rather than bool
#     zernike_indexes = centrosome.zernike.get_zernike_indexes( zernike_degree + 1)


#     z = zernike_polynomial(labels, zernike_indexes)

#     z_cm = zernike_metrics(cellmask_proj, z)
#     z_org = zernike_metrics(org_proj, z)
#     z_nuc = zernike_metrics(nucleus_proj, z)


#     # nm_labels = [f"{n}_{m}" for (n, m) in (zernike_indexes)
#     stats_tab = pd.DataFrame({'organelle':organelle_name,
#                                 'mask':'cell',
#                                 'zernike_n':[zernike_indexes[:,0].tolist()],
#                                 'zernike_m':[zernike_indexes[:,1].tolist()],
#                                 'zernike_cm_mag':[z_cm[0].tolist()],
#                                 'zernike_cm_phs':[z_cm[1].tolist()],   
#                                 'zernike_obj_mag':[z_org[0].tolist()],
#                                 'zernike_obj_phs':[z_org[1].tolist()],
#                                 'zernike_nuc_mag':[z_nuc[0].tolist()],
#                                 'zernike_nuc_phs':[z_nuc[1].tolist()]})

#     return stats_tab


# def get_XY_distribution(        
#         mask: np.ndarray,
#         centering_obj: np.ndarray,
#         obj:np.ndarray,
#         obj_name: str,
#         scale: Union[tuple, None]=None,
#         num_bins: Union[int, None] = None,
#         center_on: bool = False,
#         keep_center_as_bin: bool = True,
#         zernike_degrees: Union[int,None] = None):

#     """
#     Params
#     ----------
#     mask_obj: np.ndarray,
#         a binary 3D (ZYX) np.ndarray of the area that will be measured from
#     centering_obj: np.ndarray
#         a binary 3D (ZYX) np.ndarray of the object that will be used as the center of the concentric rins ("bins")
#     obj: np.ndarray
#         a 3D (ZYX) np.ndarray image of what will be measured within the masked area
#     obj_name: str
#         the name or nickname for the obj being measured; this will appear as a column in the output datasheet
#     scale: Union[tuple, None]=None
#         a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)
#     num_bins: Union[int,None] = None
#         the number of concentric rings to draw between the centering object and edge of the mask; None will result in 5 bins
#     center_on: bool = False
#         True = distribute the bins from the center of the centering object
#         False = distribute the bins from the edge of the centering object
#     keep_center_as_bin: bool = True
#         True = include the centering object area when creating the bins
#         False = do not include the centering object area when creating the bins
#     zernike_degrees: Union[int,None] = None
#         the number of zernike degrees to include for the zernike shape descriptors; if None, the zernike measurements will not 
#         be included in the output


#     Returns
#     -----------
#     XY_metrics:
#         a pandas Dataframe of bin, wedge, and zernike measurements
#     dist_bin_mask:
#         an np.ndarray mask of the concentric ring bins
#     dist_wedge_mask 
#         an np.ndarray mask of the 8 radial wedges

#     """

#     mask_proj = create_masked_sum_projection(mask)
#     center_proj = create_masked_sum_projection(centering_obj,mask.astype(bool))
#     obj_proj = create_masked_sum_projection(obj,mask.astype(bool))
 

#     XY_metrics, dist_bin_mask, dist_wedge_mask = get_concentric_distribution(mask_proj=mask_proj, 
#                                                         centering_proj=center_proj, 
#                                                         obj_proj=obj_proj, 
#                                                         obj_name=obj_name, 
#                                                         scale=scale,
#                                                         bin_count=num_bins, 
#                                                         center_on=center_on,
#                                                         keep_center_as_bin=keep_center_as_bin)
    
#     if zernike_degrees is not None:
#         zernike_metrics = get_zernike_metrics(cellmask_proj=mask_proj, 
#                                             org_proj=obj_proj,
#                                             organelle_name=obj_name, 
#                                             nucleus_proj=center_proj, 
#                                             zernike_degree=zernike_degrees)
        
#         XY_metrics = pd.merge(XY_metrics, zernike_metrics, on="object")

#     return XY_metrics, dist_bin_mask, dist_wedge_mask

# def get_XY_distribution(
#         cellmask_proj: np.ndarray,
#         nucleus_proj: np.ndarray,
#         org_proj: np.ndarray,
#         org_name: str,
#         bin_count: Union[int, None] = 5,
#         center_obj_as_bin: bool = True,
#         bins_from_center:bool = False
#     ):
#     """
#     Based on CellProfiler's measureobjectintensitydistribution. Measure the distribution of segmented objects within a masked area. 
#     In our case, we will usually utilize this function to measure the amount of an organelle within the cell.
#     Radial bins are created out from a center point, usually the nucleus edge.

#     Parameters
#     ------------
#     cellmask_proj: np.ndarray
#         a sum projection of the segmented cell area where the "intensity" value of each pixel is equal to the number of z slices where the binary cell mask is True
#     nucleus_proj: np.ndarray
#         a sum projection of the segmented nucleus area where the "intensity" value of each pixel is equal to the number of z slices where the binary nucleus mask is True
#     org_proj: np.ndarray,
#         a sum projection of the segmented organelle area where the "intensity" value of each pixel is equal to the number of z slices where the binary organelle mask is True
#     org_name: str,
#         the name or nickname of your organelle; used for labeling columns in the dataframe
#     bin_count: Union[int, None] = 5,
#         the number of bins to create within the cell mask
#     center_obj_as_bin: bool = True,
#         True = include the centering object area when creating the bins
#         False = do not include the centering object area when creating the bins
#     bins_from_center:bool = False
#         True = distribute the bins from the center of the centering object
#     masked


#     Returns
#     -------------
#     returns one statistics table (pd.DataFrame) + bin_array (np.ndarray) image
#     """

#     # other parameters that will stay constant
#     nobjects = 1

#     # create binary arrays
#     center_objects = nucleus_proj>0 
#     cellmask = (cellmask_proj>0).astype(np.uint16)


#     ################   ################
#     ## define masks for computing distances
#     ################   ################
#     normalized_distance, good_mask, i_center, j_center = get_normalized_distance_and_mask(cellmask, center_objects, bins_from_center, center_obj_as_bin)
    
#     if normalized_distance is None:
#         print('WTF!!  normalized_distance returned wrong')

#     ################   ################
#     ## get histograms
#     ################   ################
#     ngood_pixels = np.sum(good_mask)
#     good_labels = cellmask[good_mask]

#     # protect against None normaized_distances
#     bin_array = (normalized_distance * bin_count).astype(int)
#     bin_array[bin_array > bin_count] = bin_count # shouldn't do anything

#     #                 (    i          ,         j              )
#     labels_and_bins = (good_labels - 1, bin_array[good_mask])

#     #                coo_matrix( (             data,             (i, j)    ), shape=                      )
#     histogram_cmsk = coo_matrix( (cellmask_proj[good_mask], labels_and_bins), shape=(nobjects, bin_count) ).toarray()
#     histogram_org = coo_matrix(  (org_proj[good_mask],      labels_and_bins), shape=(nobjects, bin_count) ).toarray()

#     bin_array = (normalized_distance * bin_count).astype(int)

#     sum_by_object_cmsk = np.sum(histogram_cmsk, 1) # flattened cellmask voxel count
#     sum_by_object_org = np.sum(histogram_org, 1)  # organelle voxel count


#     # DEPRICATE: since we are NOT computing object_i by object_i (individual organelle labels)
#     # sum_by_object_per_bin = np.dstack([sum_by_object] * (bin_count + 1))[0]
#     # fraction_at_distance = histogram / sum_by_object_per_bin

#     # number of bins.
#     number_at_distance = coo_matrix(( np.ones(ngood_pixels), labels_and_bins), (nobjects, bin_count)).toarray()

#     # sicne we aren't breaking objects apart this is just ngood_pixels

#     sum_by_object = np.sum(number_at_distance, 1)

#     sum_by_object_per_bin = np.dstack([sum_by_object] * (bin_count))[0]
#     fraction_at_bin = number_at_distance / sum_by_object_per_bin # sums to 1.0

#     # object_mask = number_at_distance > 0
#     # DEPRICATE:# not doing over multiple objects so don't need object mask.. or fractionals
#     # mean_pixel_fraction = fraction_at_distance / ( fraction_at_bin + np.finfo(float).eps )
#     # masked_fraction_at_distance = np.ma.masked_array( fraction_at_distance, ~object_mask )
#     # masked_mean_pixel_fraction = np.ma.masked_array(mean_pixel_fraction, ~object_mask)

#     ################   ################
#     ## collect Anisotropy calculation.  + summarize
#     ################   ################
#     # Split each cell into eight wedges, then compute coefficient of variation of the wedges' mean intensities
#     # in each ring. Compute each pixel's delta from the center object's centroid
#     i, j = np.mgrid[0 : cellmask.shape[0], 0 : cellmask.shape[1]]
#     imask = i[good_mask] > i_center[good_mask]
#     jmask = j[good_mask] > j_center[good_mask]
#     absmask = abs(i[good_mask] - i_center[good_mask]) > abs(j[good_mask] - j_center[good_mask])
#     radial_index = (imask.astype(int) + jmask.astype(int) * 2 + absmask.astype(int) * 4)

#     # return radial_index, labels, good_mask, bin_indexes
#     stat_names =[]
#     cv_cmsk = []
#     cv_obj = []

#     # collect the numbers from each "bin"
#     for bin in range(bin_count):
#         bin_mask = good_mask & (bin_array == bin)
#         bin_pixels = np.sum(bin_mask)

#         bin_labels = cellmask[bin_mask]

#         bin_radial_index = radial_index[bin_array[good_mask] == bin]
#         labels_and_radii = (bin_labels - 1, bin_radial_index)
#         pixel_count = coo_matrix( (np.ones(bin_pixels), labels_and_radii), (nobjects, 8) ).toarray()

#         radial_counts_cmsk = coo_matrix( (cellmask_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray()
#         radial_counts = coo_matrix( (org_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray()
#         # radial_values = coo_matrix( (img_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray()

#         # we might need the masked arrays for some organelles... but I think not. keeping for now
#         mask = pixel_count == 0

#         radial_means_cmsk = np.ma.masked_array(radial_counts_cmsk / pixel_count, mask)
#         radial_cv_cmsk = np.std(radial_means_cmsk, 1) / np.mean(radial_means_cmsk, 1)
#         radial_cv_cmsk[np.sum(~mask, 1) == 0] = 0
#         radial_cv_cmsk.mask = np.sum(~mask, 1) == 0


#         radial_means_obj = np.ma.masked_array(radial_counts / pixel_count, mask)
#         radial_cv_obj = np.std(radial_means_obj, 1) / np.mean(radial_means_obj, 1)
#         radial_cv_obj[np.sum(~mask, 1) == 0] = 0
#         radial_cv_obj.mask = np.sum(~mask, 1) == 0

#         bin_name = str(bin + 1) if bin > 0 else "1"

#         stat_names.append(bin_name)
#         cv_cmsk.append(float(np.mean(radial_cv_cmsk)))  #convert to float to make importing from csv more straightforward
#         cv_obj.append(float(np.mean(radial_cv_obj)))
    
#     stats_dict={'organelle': org_name,
#                 'mask': 'cell',
#                 'radial_n_bins': bin_count,
#                 'radial_bins': [stat_names],
#                 'radial_cm_vox_cnt': [histogram_cmsk.squeeze().tolist()],
#                 'radial_org_vox_cnt': [histogram_org.squeeze().tolist()],
#                 # 'radial_org_intensity': [histogram_img.squeeze().tolist()],
#                 'radial_n_pix': [number_at_distance.squeeze().tolist()],
#                 'radial_cm_cv':[cv_cmsk],
#                 'radial_org_cv':[cv_obj]}

#     # stats_tab = pd.DataFrame(statistics,columns=col_names)
#     stats_tab = pd.DataFrame(stats_dict)  
#     return stats_tab, bin_array

# def get_XY_distribution(
#         cellmask_proj: np.ndarray,
#         org_proj: np.ndarray,
#         img_proj: np.ndarray,
#         org_name: str,
#         nucleus_proj: np.ndarray,
#         n_bins: int = 5,
#         from_edges: bool = True,
#     ):
#     """Perform the radial measurements on the image set

#     Parameters
#     ------------
#     cellmask_proj: np.ndarray,
#     org_proj: np.ndarray,
#     img_proj: np.ndarray,
#     org_name: str,
#     nucleus_proj: Union[np.ndarray, None],
#     n_bins: int = 5,
#     from_edges: bool = True,

#     masked

#     # params
#     #   n_bins .e.g. 6
#     #   normalizer - cellmask_voxels, organelle_voxels, cellmask_and_organelle_voxels
#     #   from_edges = True


#     Returns
#     -------------
#     returns one statistics table + bin_indexes image array
#     """

#     # other params
#     bin_count = n_bins if n_bins is not None else 5
#     nobjects = 1
#     scale_bins = True 
#     keep_nuc_bins = True # this toggles whether to count things inside the nuclei mask.  
#     center_on_nuc = False # choosing the edge of the nuclei or the center as the center to propogate from

#     center_objects = nucleus_proj>0 

#     # labels = label(cellmask_proj>0) #extent as 0,1 rather than bool    
#     labels = (cellmask_proj>0).astype(np.uint16)
#     # labels = np.zeros_like(cellmask_proj)
#     # labels[labels>0]=1

#     ################   ################
#     ## define masks for computing distances
#     ################   ################
#     normalized_distance, good_mask, i_center, j_center = get_normalized_distance_and_mask(labels, center_objects, center_on_nuc, keep_nuc_bins)
    
#     if normalized_distance is None:
#         print('WTF!!  normailzed_distance returned wrong')

#     ################   ################
#     ## get histograms
#     ################   ################
#     ngood_pixels = np.sum(good_mask)
#     good_labels = labels[good_mask]

#     # protect against None normaized_distances
#     bin_indexes = (normalized_distance * bin_count).astype(int)
#     bin_indexes[bin_indexes > bin_count] = bin_count # shouldn't do anything

#     #                 (    i          ,         j              )
#     labels_and_bins = (good_labels - 1, bin_indexes[good_mask])

#     #                coo_matrix( (             data,             (i, j)    ), shape=                      )
#     histogram_cmsk = coo_matrix( (cellmask_proj[good_mask], labels_and_bins), shape=(nobjects, bin_count) ).toarray()
#     histogram_org = coo_matrix(  (org_proj[good_mask],      labels_and_bins), shape=(nobjects, bin_count) ).toarray()
#     histogram_img = coo_matrix(  (img_proj[good_mask],      labels_and_bins), shape=(nobjects, bin_count) ).toarray()

#     bin_indexes = (normalized_distance * bin_count).astype(int)

#     sum_by_object_cmsk = np.sum(histogram_cmsk, 1) # flattened cellmask voxel count
#     sum_by_object_org = np.sum(histogram_org, 1)  # organelle voxel count
#     sum_by_object_img = np.sum(histogram_img, 1)  # image intensity projection

#     # DEPRICATE: since we are NOT computing object_i by object_i (individual organelle labels)
#     # sum_by_object_per_bin = np.dstack([sum_by_object] * (bin_count + 1))[0]
#     # fraction_at_distance = histogram / sum_by_object_per_bin

#     # number of bins.
#     number_at_distance = coo_matrix(( np.ones(ngood_pixels), labels_and_bins), (nobjects, bin_count)).toarray()

#     # sicne we aren't breaking objects apart this is just ngood_pixels

#     sum_by_object = np.sum(number_at_distance, 1)

#     sum_by_object_per_bin = np.dstack([sum_by_object] * (bin_count))[0]
#     fraction_at_bin = number_at_distance / sum_by_object_per_bin # sums to 1.0

#     # object_mask = number_at_distance > 0
#     # DEPRICATE:# not doing over multiple objects so don't need object mask.. or fractionals
#     # mean_pixel_fraction = fraction_at_distance / ( fraction_at_bin + np.finfo(float).eps )
#     # masked_fraction_at_distance = np.ma.masked_array( fraction_at_distance, ~object_mask )
#     # masked_mean_pixel_fraction = np.ma.masked_array(mean_pixel_fraction, ~object_mask)

#     ################   ################
#     ## collect Anisotropy calculation.  + summarize
#     ################   ################
#     # Split each cell into eight wedges, then compute coefficient of variation of the wedges' mean intensities
#     # in each ring. Compute each pixel's delta from the center object's centroid
#     i, j = np.mgrid[0 : labels.shape[0], 0 : labels.shape[1]]
#     imask = i[good_mask] > i_center[good_mask]
#     jmask = j[good_mask] > j_center[good_mask]
#     absmask = abs(i[good_mask] - i_center[good_mask]) > abs(
#         j[good_mask] - j_center[good_mask]
#     )
#     radial_index = (
#         imask.astype(int) + jmask.astype(int) * 2 + absmask.astype(int) * 4
#     )

#     # return radial_index, labels, good_mask, bin_indexes
#     statistics = []
#     stat_names =[]
#     cv_cmsk = []
#     cv_obj = []
#     cv_img = []
#     # collect the numbers from each "bin"
#     for bin in range(bin_count):
#         bin_mask = good_mask & (bin_indexes == bin)
#         bin_pixels = np.sum(bin_mask)

#         bin_labels = labels[bin_mask]

#         bin_radial_index = radial_index[bin_indexes[good_mask] == bin]
#         labels_and_radii = (bin_labels - 1, bin_radial_index)
#         pixel_count = coo_matrix( (np.ones(bin_pixels), labels_and_radii), (nobjects, 8) ).toarray()

#         radial_counts_cmsk = coo_matrix( (cellmask_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray()
#         radial_counts = coo_matrix( (org_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray()
#         radial_values = coo_matrix( (img_proj[bin_mask], labels_and_radii), (nobjects, 8) ).toarray()

#         # we might need the masked arrays for some organelles... but I think not. keeping for now
#         mask = pixel_count == 0

#         radial_means_cmsk = np.ma.masked_array(radial_counts_cmsk / pixel_count, mask)
#         radial_cv_cmsk = np.std(radial_means_cmsk, 1) / np.mean(radial_means_cmsk, 1)
#         radial_cv_cmsk[np.sum(~mask, 1) == 0] = 0
#         radial_cv_cmsk.mask = np.sum(~mask, 1) == 0


#         radial_means_obj = np.ma.masked_array(radial_counts / pixel_count, mask)
#         radial_cv_obj = np.std(radial_means_obj, 1) / np.mean(radial_means_obj, 1)
#         radial_cv_obj[np.sum(~mask, 1) == 0] = 0
#         radial_cv_obj.mask = np.sum(~mask, 1) == 0

#         radial_means_img = np.ma.masked_array(radial_values / pixel_count, mask)
#         radial_cv_img = np.std(radial_means_img, 1) / np.mean(radial_means_img, 1)
#         radial_cv_img[np.sum(~mask, 1) == 0] = 0
#         radial_cv_img.mask = np.sum(~mask, 1) == 0

#         bin_name = str(bin) if bin > 0 else "Ctr"

#         # # there's gotta be a better way to collect this stuff together... pandas?
#         # statistics += [
#         #     (   bin_name,
#         #         # np.mean(number_at_distance[:, bin]), 
#         #         # np.mean(histogram_cmsk[:, bin]), 
#         #         # np.mean(histogram_org[:, bin]), 
#         #         # np.mean(histogram_img[:, bin]), 
#         #         np.mean(radial_cv_cmsk) ,
#         #         np.mean(radial_cv_obj) ,
#         #         np.mean(radial_cv_img) )
#         # ]
#         stat_names.append(bin_name)
#         cv_cmsk.append(float(np.mean(radial_cv_cmsk)))  #convert to float to make importing from csv more straightforward
#         cv_obj.append(float(np.mean(radial_cv_obj)))
#         cv_img.append(float(np.mean(radial_cv_obj)))

#     # TODO: fix this grooooos hack
#     # col_names=['organelle','mask','bin','n_bins','n_pix','cm_vox_cnt','org_vox_cnt','org_intensity','cm_radial_cv','org_radial_cv','img_radial_cv']
#     # stats_dict={'organelle': org_name,
#     #             'mask': 'cell',
#     #             'radial_n_bins': bin_count,
#     #             'radial_bins': [[s[0] for s in statistics]],
#     #             'radial_cm_vox_cnt': [histogram_cmsk.squeeze().tolist()],
#     #             'radial_org_vox_cnt': [histogram_org.squeeze().tolist()],
#     #             'radial_org_intensity': [histogram_img.squeeze().tolist()],
#     #             'radial_n_pix': [number_at_distance.squeeze().tolist()],
#     #             'radial_cm_cv':[[s[1] for s in statistics]],
#     #             'radial_org_cv':[[s[2] for s in statistics]],
#     #             'radial_img_cv':[[s[3] for s in statistics]],
#     #             }
    
#     stats_dict={'organelle': org_name,
#                 'mask': 'cell',
#                 'radial_n_bins': bin_count,
#                 'radial_bins': [stat_names],
#                 'radial_cm_vox_cnt': [histogram_cmsk.squeeze().tolist()],
#                 'radial_org_vox_cnt': [histogram_org.squeeze().tolist()],
#                 'radial_org_intensity': [histogram_img.squeeze().tolist()],
#                 'radial_n_pix': [number_at_distance.squeeze().tolist()],
#                 'radial_cm_cv':[cv_cmsk],
#                 'radial_org_cv':[cv_obj],
#                 'radial_img_cv':[cv_img],
#                 }

#     # stats_tab = pd.DataFrame(statistics,columns=col_names)
#     stats_tab = pd.DataFrame(stats_dict)  
#     return stats_tab, bin_indexes

# quantify the distribution of one or more organelles from one cell
# USED #
def get_distribution_metrics(source_file_path: str,
                        list_obj_names: List[str],
                        list_obj_segs: List[np.ndarray],
                        list_region_names: Union[List[str], None]=None,
                        list_region_segs: Union[List[np.ndarray], None]=None,
                        mask_name: Union[str, None]=None,
                        scale: Union[tuple,None] = None,
                        centering_obj: Union[str, None]=None,
                        num_bins: Union[int, None]=5,
                        center_on: Union[bool, None]=False,
                        keep_center_as_bin: Union[bool, None]=True,
                        zernike_degrees: Union[int, None]=9) -> pd.DataFrame:
    """
    Measure the spatial distribution of multiple organelles from a single cell in respect to the center

    Parameters:
    ----------
    source_file: str
        Path to the source image file. This will be used as part of the metadata information in the output table. 
    list_obj_names: List[str]
        List of organelle names. These names should match the suffix on the segmentation image files.
    list_obj_segs: List[np.ndarray]
        List of 3D organelle segmentation arrays matching the order included in list_obj_names.
    list_region_names: Union[List[str], None]
        List of segmented region/mask names. These names should match the suffix on the segmentation image files.
        This should include:
            - a mask segmentation, such as the cell mask, for masking during distribution analysis; else, the entire image will be 
            quantified. Only one object per mask image will be analyzed. If there are more than one included, they will be combined 
            prior to analysis and the entire region will be quantified. If no mask is provided, the entire image will be quantified.
            - a centering object, such as the nucleus; else the center of the mask region will be used as the XY distribution centering 
            point.
    list_region_segs: Union[List[np.ndarray], None]
        List of 3D region segmentation arrays matching the order specified in list_region_names. Specify None if no regions are provided.
    mask_name: Union[str, None]
        Name of the region to use as the mask for analysis; if not specified, the entire image will be quantified.
    scale: Union[tuple,None] = None
        a tuple that contains the real world dimensions for each dimension in the image (Z, Y, X)
    centering_obj : Union[str, None], default=None
        Name of the region to use for centering during distribution analysis.
        This region should be included in the list_region_names and list_region_segs variables.
        If not specified, the center of the mask, or entire image if no mask was specified, will be used as the centering object.
    num_bins : Union[int, None], default=5
        Number of radial bins to create in the XY distribution analysis.
    center_on : Union[bool, None], default=True
        Whether to start creation of the XY bins from the center (True) or the edge (False) of the centering object.
    keep_center_as_bin : Union[bool, None], default=True
        Whether to keep the centering object as the first XY bin. 
    zernike_degrees : Union[int, None], default=9
        Zernike polynomial degree for circular shape/pattern analysis in the XY distribution analysis.
        If None and include_dist=True, no Zernike features will be calculated.

    Returns
    -------
    final_dist_tab : pd.DataFrame
        Dataframe for XY and Z distribution metrics for organelle in one image

    """
    # Validate inputs
    if not list_obj_names:
        raise ValueError("list_obj_names cannot be empty")
    if len(list_obj_names) != len(list_obj_segs):
        raise ValueError(f"Mismatch: {len(list_obj_names)} items in list_obj_names but {len(list_obj_segs)} items in list_obj_segs")
    if list_region_names and list_region_segs and len(list_region_names) != len(list_region_segs):
        raise ValueError(f"Mismatch: {len(list_region_names)} items in list_region_names but {len(list_region_segs)} items in list_region_segs")
    
    if isinstance(source_file_path, str): source_file_path = Path(source_file_path)
    print(f"Quantifying organelle interactions from {source_file_path.name}")

    # specify the mask image to use during quantification based on the mask_name provided
    if list_region_names is None or list_region_segs is None:
        print("No regions provided. No mask or centering object will be applied before analysis.")
        mask = None
        centering_obj = None
    elif mask_name is None or mask_name not in list_region_names:
        if mask_name is not None:
            raise ValueError(f"Mask '{mask_name}' not found. No mask will be applied before analysis.")
        mask = None
        mask_name = None
    else:
        mask = list_region_segs[list_region_names.index(mask_name)]

    mask_name = "whole_image" if mask_name is None else mask_name

    # specify the centering image to use during quantification based on the centering object name provided
    if centering_obj == None:
        print("No centering object provided. Using center of mask or entire image for distribution centering.")
        centering_img = None
    elif centering_obj not in list_region_names:
        raise ValueError(f"Centering object '{centering_obj}' not found in region names: {list_region_names}")
    else:
        centering_img = list_region_segs[list_region_names.index(centering_obj)]
    
    # empty list to collect the distribution data for each organelle
    dist_tabs = []
    XY_bins_imgs = []
    XY_wedges_imgs = []

    # loop through the list of organelles and run the get_XY_distribution and get_Z_distribution function
    for j, target in enumerate(list_obj_names):    
        # select segmentation and if ER, ensure it is only one object
        if target == 'ER':
            org_obj = (list_obj_segs[j] > 0).astype(np.uint16)
        else:
            org_obj = list_obj_segs[j]

        # run get_XY_distribution function to output a table of distribution measurements in respect to a specified object in the XY
        XY_distribution, XY_bins, XY_wedges = get_XY_distribution(mask=mask,
                                                                    mask_name = mask_name,
                                                                    centering_obj=centering_img,
                                                                    obj=org_obj,
                                                                    obj_name=target,
                                                                    scale=scale,
                                                                    num_bins=num_bins,
                                                                    center_on=center_on,
                                                                    keep_center_as_bin=keep_center_as_bin,
                                                                    zernike_degrees=zernike_degrees # set to None if you wish to skip quantification of zernike features
                                                                    )
        
        # if XY_bins_imgs list is empty append, if not skip
        if not XY_bins_imgs and not XY_wedges_imgs:
            XY_bins_imgs.append(XY_bins)
            XY_wedges_imgs.append(XY_wedges)
            
        # run get_Z_distribution function to output a table of distribution measurements in respect to a specified object in the Z
        Z_distribution = get_Z_distribution(mask=mask,
                                                mask_name = mask_name, 
                                                obj=org_obj,
                                                obj_name=target,
                                                center_obj=centering_img,
                                                scale=scale)

        # add table to list above
        dist_tab = pd.merge(XY_distribution, Z_distribution)
        dist_tabs.append(dist_tab)

    # combine the lists for each organelle into one table
    final_dist_tab = pd.concat(dist_tabs, ignore_index=True)

    # add a new column to list the name of the image these data are derived from 
    final_dist_tab.insert(loc=0,column='image_name',value=source_file_path.stem)

    return final_dist_tab

# batch process distribution quantification for multiple cells from a single experiment
def batch_process_distribution_quant(dataset_name: str,
                             raw_path: Union[Path,str], 
                             seg_path: Union[Path,str],
                             quant_path: Union[Path, str], 
                             raw_file_type: str,
                             organelle_names: List[str],
                             region_names: Union[List[str], None]=None,
                             mask_name: Union[str, None]=None,
                             use_scale:bool=True,
                             seg_suffix:Union[str, None]=None,
                             centering_obj: Union[str, None]=None,
                             num_bins: Union[int, None]=5,
                             center_on: Union[bool, None]=False,
                             keep_center_as_bin: Union[bool, None]=True,
                             zernike_degrees: Union[int, None]=9):
    """  
    batch process distribution quantification; this function is currently optimized to process images from one file folder per image type (e.g., raw, segmentation)
    the output csv files are saved to the indicated quant_path folder

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
    organelle_names: List[str]
        List of organelle names to analyze. These names should match the suffix on the organelle segmentation files
    region_names: Union[List[str], None]=None
        List of region names to analyze. Usually ['cell', 'nuc'] for cell mask and nucleus.
        If no regions are to be included, specify None here.
    mask_name: Union[str, None]=None
        Name of the region to use for segmentation (if any). This name should be included in the regions_name variable.
        If None, the entire image will be quantified.
    use_scale: bool=True
        Whether to apply scaling to the quantitative data; scaled data will be in real world units (e.g., microns) rather than pixels/voxels
    seg_suffix:Union[str, None]=None
        Any additional text that is included in the segmentation tiff files between the file stem and the segmentation suffix, not including the initial "-"
    centering_obj : str or None, default=None
        Name of the region to use for centering distribution analysis
        This region should be included in the list_region_names and list_region_segs variables
        If not specified, the center of the mask, or entire image if no mask was specified, will be used as the centering object
    num_bins : int or None, default=5
        Number of radial bins to create in the XY distribution analysis
    center_on : bool or None, default=True
        Whether to start creation of the XY bins from the center of the centering object (True) or edge (False)
    keep_center_as_bin : bool or None
        Whether to keep centering object as the first XY bin
    zernike_degrees : int or None, default=9
        Zernike polynomial degree for shape analysis in the XY distribution analysis
        If None and include_dist=True, no Zernike features will be calculated
    
    Returns:
    --------
    None
        Saves output files to the specified quantification path
    """

    start = time.time()
    count = 0

    # create path objects if inputs are strings
    if isinstance(raw_path, str): raw_path = Path(raw_path)
    if isinstance(seg_path, str): seg_path = Path(seg_path)
    if isinstance(quant_path, str): quant_path = Path(quant_path)
    
    # create directory if it doesn't exist
    if not Path.exists(quant_path):
        Path.mkdir(quant_path)
        print(f"Output file path not found. Making {quant_path}.")
    
    # specify the columns that will be checked per table
    unique_keys = ['dataset', 'image_name']

    # check if any existing data is present in outfiles
    dist_path = quant_path / f"{dataset_name}_distribution_metrics.csv"
    existing_dist_keys = load_existing_keys_csv(dist_path, unique_keys)

    # list of organelle segmentation and masks files to collect from each image
    segs_to_collect = organelle_names + region_names if region_names is not None else organelle_names

    # reading list of files from the raw path
    img_file_list = list_image_files(raw_path, raw_file_type)
    len_file_list = len(img_file_list)

    # loop through list of images and quantify distribution metrics if data for that image do not already exist;
    for img_f in img_file_list:
        img_start = time.time()
        count = count + 1
        # skip files that have already been processed
        if (dataset_name, img_f.stem) in existing_dist_keys:
            print(f"Skipping {img_f.name} as it is already listed in the output file(s).")
            continue

        # process analysis for this cells
        else:
            filez = find_segmentation_tiff_files(img_f, segs_to_collect, seg_path, seg_suffix)

            # read in raw file and metadata
            _img_data, meta_dict = read_czi_image(filez["raw"])

            # store organelle images as list
            organelles = [read_tiff_image(filez[org]) for org in organelle_names]

            # load regions as a list based on order in list (should match order in "masks" file)
            regions = [read_tiff_image(filez[r]) for r in region_names] if region_names is not None else None

            # define the scale
            if use_scale is True:
                scale_tup = meta_dict['scale']
            else:
                scale_tup = None

            dist_tab = get_distribution_metrics(source_file_path=img_f,
                                            list_obj_names=organelle_names,
                                            list_obj_segs=organelles, 
                                            list_region_names=region_names,
                                            list_region_segs=regions, 
                                            mask_name=mask_name,
                                            scale=scale_tup,
                                            centering_obj=centering_obj,
                                            num_bins=num_bins,
                                            center_on=center_on,
                                            keep_center_as_bin=keep_center_as_bin,
                                            zernike_degrees=zernike_degrees)
            
            dist_tab = dist_tab.astype(str)  # ensure all data is string to avoid dtype issues

            # save the distribution (or labels only) table data per image directly to csv
            dist_tab.insert(loc=0,column='dataset',value=dataset_name)
            append_atomic_csv(dist_path, dist_tab)
            del dist_tab  # free up memory

            end2 = time.time()
            print(f"Completed distribution quantification of {meta_dict['file_name']} in {(end2-img_start)/60} mins.")
            print(f"{count}/{len_file_list} images have been processed.")
            print(f"Time elapsed: {(end2-img_start)/60} mins")

    end = time.time()
    print(f"Distribution quantification for {count} files is COMPLETE! Files saved to '{quant_path}'.")
    print(f"It took {(end - start)/60} minutes to quantify these files.")

# summarize distribution values per organelle per cell across one or more experiments
def batch_distribution_summary_stats(out_prefix: str,
                                      csv_path_list: List[str],
                                      out_path: str,
                                      mask_name: str = "whole_image"):
    """ 
    Batch process interaction quantification summary statistics from multiple datasets.

    Parameters:
    -----------
    out_prefix: str
        The prefix used to name the output file. An "_" will be included between this prefix and the file suffix.
    csv_path_list: List[str],
        A list of path strings where .csv files to analyze are located.
    out_path: str,
        A path string where the summary data file will be output to
    mask_name: str = "whole_image"
        Name of the region to use as the mask for analysis across all datasets
    """

    # for keeping track of dataset and file numbers
    ds_count = 0
    fl_count = 0

    ###################
    # Read in the csv files and combine them into one
    ###################
    # create empty list to hold the distribution tables from different experiments
    dist_tabs = []

    # loop through all of the locations listed above and find the _distribution files; append them to the list above
    for loc in csv_path_list:

        # list all csv files in the location
        files_store = sorted(loc.glob("*.csv"))

        # find the unique datasets in this location based on the prefixes before "_distribution_metrics"
        prefixes = set(f.name.split("_distribution_metrics")[0] for f in files_store if "_distribution_metrics" in f.name)
        print(f"Found the following datasets in {loc}:", prefixes)
        for prefix in prefixes:
            ds_count += 1
            # select only the files from this dataset
            files_subset = [f for f in files_store if f.name.startswith(prefix +"_distribution_metrics")]
            for file in files_subset:
                fl_count += 1
                stem = file.stem
                if "_dist" in stem:
                    test_dist = pd.read_csv(file, index_col=0)
                    dist_tabs.append(test_dist)

    # combine the dist lists found above into one table
    dist_df = pd.concat(dist_tabs,axis=0, join='outer').reset_index()

    print(f"Found {fl_count} files from {ds_count} dataset(s) across {len(csv_path_list)} location(s).")

    # mask name checker
    mask_name = "whole_image" if mask_name is None else mask_name
    
    # extract centering object metrics
    if 'XY_center_vox_cnt_perbin' in list(dist_df.columns): # if there is a centering object
        nuc_dist_df = dist_df[["dataset", "image_name", 'scale',
                            "XY_bins", "XY_center_vox_cnt_perbin", f"XY_{mask_name}_vox_cnt_perbin", "XY_center_cv_perbin",
                            "XY_wedges", "XY_center_vox_cnt_perwedge", f"XY_{mask_name}_vox_cnt_perwedge",
                            "Z_slices", "Z_center_vox_cnt", f"Z_{mask_name}_vox_cnt"]].drop_duplicates(subset=['dataset', 'image_name'])
        nuc_dist_df.columns = nuc_dist_df.columns.str.replace('center', 'obj', regex=False)
        nuc_dist_df.insert(loc=3,column='object',value='nuc')
        nuc_dist_df.set_index(['dataset', 'image_name', 'scale', 'object'], inplace=True)

        # select relevant columns from dist dataset
        dist_df2 = dist_df[list(nuc_dist_df.reset_index().columns)]
        dist_df2.set_index(['dataset', 'image_name', 'scale', 'object'], inplace=True)

        # combine
        combo_dist_df = pd.concat([nuc_dist_df, dist_df2], axis=0)
    else: # if there is not a centering object
        dist_df.set_index(['dataset', 'image_name', 'scale', 'object'], inplace=True)
        combo_dist_df = dist_df

    # loop through each row of data and calculate histogram statistics
    hist_dfs = []
    for ind in combo_dist_df.index:
        selection = combo_dist_df.loc[[ind]].reset_index()
        bins_df = pd.DataFrame()
        wedges_df = pd.DataFrame()
        Z_df = pd.DataFrame()
        CV_df = pd.DataFrame()

        # select relevant columns into different groups
        bins_df[['bins', 'masks', 'obj']] = selection[['XY_bins', f'XY_{mask_name}_vox_cnt_perbin', 'XY_obj_vox_cnt_perbin']]
        wedges_df[['bins', 'masks', 'obj']] = selection[['XY_wedges', f'XY_{mask_name}_vox_cnt_perwedge', 'XY_obj_vox_cnt_perwedge']]
        Z_df[['bins', 'masks', 'obj']] = selection[['Z_slices', f'Z_{mask_name}_vox_cnt', 'Z_obj_vox_cnt']]
        CV_df[['XY_obj_cv_perbin']] = selection[['XY_obj_cv_perbin']]

        dfs = [selection[['dataset', 'image_name', 'scale', 'object']].reset_index()]

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
        
        # combine dataframes per group together
        combined_df = pd.concat(dfs,axis=1).drop(columns="index")
        combined_df.set_index(['dataset', 'image_name', 'scale', 'object'], inplace=True)
        hist_dfs.append(combined_df)

    # combine data from each row of data in the original table together
    dist_summary = pd.concat(hist_dfs).sort_values(by=['dataset', 'image_name', 'scale', 'object'])
    dist_summary.reset_index(inplace=True)

    # export before unstacking
    if (Path(out_path) / f"{out_prefix}_per_org_distribution_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}_per_org_distribution_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
    else:
        dist_summary.to_csv(str(out_path) + f"/{out_prefix}_per_org_distribution_summarystats.csv")
    
    # unstack and format interaction distribution summary table
    dist_summary.insert(2, "mask_name", mask_name) ## TODO: change after dist is updated to include mask_name
    dist_final = dist_summary.set_index(['dataset', 'image_name', 'mask_name', 'scale', 'object']).unstack(-1)
    dist_final.columns = ["_".join((col_name[1], col_name[0])) for col_name in dist_final.columns.to_flat_index()]
    dist_final = dist_final.reset_index()

    ###################
    # export summary sheets
    ###################
    if (Path(out_path) / f"{out_prefix}_distribution_summarystats.csv").exists():
        raise FileExistsError(f"CAUTION: {out_prefix}_distribution_summarystats.csv already exists and will not be overwritten. Move the existing file, change the `out_prefix` or `out_path` to continue without error.")
    else:
        dist_final.to_csv(str(out_path) + f"/{out_prefix}_distribution_summarystats.csv", mode='x')
        print(f"Exported distribution summary statistics (after unstacking) to {out_path}/{out_prefix}_distribution_summarystats.csv")
    print(f"Organelle distribution summary is complete.")
    return dist_final