import numpy as np
from scipy import ndimage
from skimage.morphology import opening
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed
from aicssegmentation.core.utils import size_filter
from infer_subc.core.img import label_uint16
from skimage.measure._label import label

def _highpass_filter(in_img: np.ndarray, sigma:float=0.0,
                     iterations:int=1, open:bool=False) -> np.ndarray:
    highpass = in_img.copy().astype(np.float32)
    for _ in range(iterations):
        lowpass = ndimage.gaussian_filter(highpass, sigma)
        highpass -= lowpass
        np.clip(highpass, 0, None, out=highpass)
    if open:
        highpass=opening(highpass)
    return highpass

def _otsu_size_filter(in_img: np.ndarray, thresh_adj:float=1, min_size:int=0) -> np.ndarray:
    threshold = threshold_otsu(in_img)
    ots = (in_img >= (threshold*thresh_adj))
    return size_filter(img=ots, min_size=min_size, method='3D')

def watershed_declumping(raw_img:np.ndarray, seg_img:np.ndarray, declump:bool, 
                         sigma:float, iterations:int=1, open:bool=False, 
                         thresh_adj:float=1, min_size:int=0) -> np.ndarray:
    seg_img = label_uint16(seg_img)
    if declump and iterations>=1:
        highpass = _highpass_filter(in_img=raw_img, sigma=sigma, open=open, iterations=iterations)
        ots = _otsu_size_filter(in_img=highpass, thresh_adj=thresh_adj, min_size=min_size)
        
        # Changed to match the method_notebook
        # return label(watershed(image=(np.max(raw_img)-raw_img), 
        #                                               markers=label_uint16(ots), 
        #                                               mask=seg_img,
        #                                               connectivity=np.ones((3, 3, 3), bool))).astype(np.uint16)
        return label((seg_img) + watershed(image=(np.max(raw_img)-raw_img), 
                                                      markers=label_uint16(ots), 
                                                      mask=seg_img,
                                                      connectivity=np.ones((3, 3, 3), bool))).astype(np.uint16)
    else:
        return seg_img