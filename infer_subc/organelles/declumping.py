import numpy as np
from scipy import ndimage
from skimage.morphology import opening
from skimage.filters import threshold_otsu
from skimage.segmentation import watershed
from aicssegmentation.core.utils import size_filter
from skimage.measure._label import label

def _highpass_filter(in_img: np.ndarray, sigma:float=0.0,
                     iterations:int=1, open:bool=False) -> np.ndarray:
    """
    Applies a highpass filter to the input image.

    Parameters:
    ---------

    in_img:np.ndarray:
        The input image to apply the filter to
    sigma:float
        The sigma value used for the lowpass filter
    iterations:int
        The number of times to apply the highpass filter
    open:bool
        A true/false statement of whether to apply an opening filter or not
    """
    highpass = in_img.copy().astype(np.float32)
    for _ in range(iterations):
        lowpass = ndimage.gaussian_filter(highpass, sigma)
        highpass -= lowpass
        np.clip(highpass, 0, None, out=highpass)
    if open:
        highpass=opening(highpass)
    return highpass

def _otsu_size_filter(in_img: np.ndarray, thresh_adj:float=1, min_size:int=0) -> np.ndarray:
    """
    Both thresholds the image and ensures that all objects are above specified size

    Parameters:
    ---------
    in_img:np.ndarray
        The input image to apply the filter to
    thresh_adj:float
        A scalar value to adjust the threshold value by
    min_size:int
        The minimum value for the volume of the objects
    
    Returns:
    ----------
    A thresholded np.ndarray
    """
    threshold = threshold_otsu(in_img)
    ots = (in_img >= (threshold*thresh_adj))
    return size_filter(img=ots, min_size=min_size, method='3D')

def watershed_declumping(raw_img:np.ndarray, seg_img:np.ndarray, declump:bool, 
                         sigma:float, iterations:int=1, open:bool=False, 
                         thresh_adj:float=1, min_size:int=0) -> np.ndarray:
    """
    Declumps the input organelle using a highpass and threshold to develop seeds for watershedding

    Parameters:
    ---------
    raw_img:np.ndarray
        The raw image to determine the peak points of intensity gradient
    seg_img:np.ndarray
        The segmentation image to determine the area the organelle is located
    declump:bool
        A true/false statement of whether to declump the organelle or not
    sigma:float
        The sigma value used in the gaussian blur lowpass for the highpass filter
    iterations:int
        The number of times the highpass filter is applied
    open:bool
        A true/false statement of whether to apply an opening filter in the highpass filter
    thresh_adj:float
        A scalar for the otsu threshold
    min_size:int
        A number used to ensure the minimum size of the "seeds" for watershedding
    
    Returns:
    ---------
    A labeled np.ndarray of individual organelles
    """
    seg_img = label(seg_img)
    if declump and iterations>=1:
        highpass = _highpass_filter(in_img=raw_img, sigma=sigma, open=open, iterations=iterations)
        ots = _otsu_size_filter(in_img=highpass, thresh_adj=thresh_adj, min_size=min_size)
        
        return label((seg_img) + watershed(image=(np.max(raw_img)-raw_img), 
                                                      markers=label(ots), 
                                                      mask=seg_img,
                                                      connectivity=np.ones((3, 3, 3), bool)))
    else:
        return seg_img