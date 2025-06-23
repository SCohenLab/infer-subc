import numpy as np
from infer_subc.core.img import scale_and_smooth
from skimage.morphology import opening
from skimage.filters import threshold_otsu
from skimage.measure import label
from skimage.segmentation import watershed
from aicssegmentation.core.utils import size_filter

def _highpass_filter(in_img: np.ndarray, sigma:float=0.0, median_size:int=0, 
                     iterations:int=1, open:bool=False) -> np.ndarray:
    highpass = in_img.copy()
    for _ in range(iterations):
        lowpass = scale_and_smooth(highpass,
                                   median_size = median_size, 
                                   gauss_sigma = sigma)
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
                         sigma:float, median_size:int, iterations:int=1, open:bool=False, 
                         thresh_adj:float=1, min_size:int=0) -> np.ndarray:
    if declump and iterations>=1:
        highpass = _highpass_filter(in_img=raw_img, sigma=sigma, median_size=median_size, 
                                    open=open, iterations=iterations)
        ots = _otsu_size_filter(in_img=highpass, thresh_adj=thresh_adj, min_size=min_size)
        return label((seg_img[seg_img>0]) + watershed(image=(np.max(raw_img)-raw_img), 
                                                      markers=label(ots), 
                                                      mask=seg_img[seg_img>0],
                                                      connectivity=np.ones((3, 3, 3), bool)))
    else:
        return seg_img