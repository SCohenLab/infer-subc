# **0.1.masks (Default “1”)**

The **0.1.masks workflow** produces instance segmentations of the **nucleus, cellmask (cell area)**, and **cytoplasm** via fluorescence confocal microscopy images. Using the workflow steps, you can adjust the parameters described below to improve the quality of your segmentations. **0.1.masks** is the default *“masks”* workflow for cells with a nucleus marker. If you are unsure that this workflow applies to your images, reference the imaging requirements listed below.

| **Imaging Requirements**               |   "1"    |    A     |    B     | 
| :------------------------------------- | :------: | :------: | :------: |
| Nuclei Marker                          |   ✔     |   ✘      |   ✘      |
| Cell Membrane Marker                   |   ✘     |   ✘      |   ✘      |
| Cytoplasmic Organelles                 |   ✔     |   ✔      |   ✔      |
| Cells in FOV                           | Multiple | Single   | Multiple  |

⚠️ Make sure that the first character of the microscopy image layer name is not numeric


## EXTRACTION

### **0.1.1 | Segment nuclei *(from label)***

In this step, a segmentation is generated for the nuclei which will be used to make the cellmask and cytoplasm mask in later steps. First, the intensities of the image are rescaled so that the minimum is 0 and the highest intensity is 1. A median and/or Gaussian filter is then applied to the nucleus channel to reduce noise. Next, a threshold is calculated and applied to the smoothed image using a logarithmic transformation and Li's minimum cross-entropy method. To conclude the nuclei segmentation, holes and objects are removed for final touches.

Parameters:

-	`nuc_chan`: the index of the nuclei channel *(Please note that the first channel is indexed as 0)*
-	`median_size`: the radius of the median filter; a large radius produces a more blurred image.
-	`gauss_sigma`:  the size of the Gaussian kernel; a larger sigma value produces a more blurred image.
-	`thresh_factor`: The ratio to be applied to the initial threshold value
-	`thresh_min`: The minimum value that can be used as a threshold
-	`thresh_max`: The maximum value that can be used as a threshold
-	`hole_min`: the minimum hole width to be filled in the segmentation post-thresholding
-	`hole_max`: the maximum hole width to be filled in the segmentation post-thresholding
-	`min_size`: the minimum width of the objects to remain in the segmentation post-thresholding
-	`fill_filter_method`: determines if both the hole filling and size filtering methods are applied in 3D or on a 2D level *(slice by slice)*

❗ If the fill_filter_method drop-down menu is set to *“slice_by_slice”* the min_size, hole_min, and hole_max values are squared, if set to *“3D”* the values are cubed.

### **0.1.2 | Create composite image**

This step combines the intensities of the organelle channels from the input image to create a composite image of summed intensity values.

Parameters:

-	`weight_ch`: the scale factor applied to the channel of the corresponding index before summation. Larger weight_ch values lead to increased influence of the organelle’s intensities in the composite image
-	`rescale`: whether the summed intensities of the composite image are normalized to a 0 to 1 scale

NOTE: To prevent a specific organelle’s intensities from affecting the composite image, set the weight_ch parameter for that organelle’s index to 0.


## PREPROCESSING

### **0.1.3 | Rescale and smooth image**

This step first rescales the composite image’s intensity values from 0 to 1. Next a median and/or Gaussian filter (for noise reduction) to the composite image.

Parameters:

-	`median_size`: the radius of the median filter; a large radius produces a more blurred image.
-	`gauss_sigma`:  the size of the Gaussian kernel; a larger sigma value produces a more blurred image.

❗ If you do not want to apply smoothing, set the parameter to 0.

### **0.1.4 | Log transform + Scharr edge detection**

This step combines two effects to emphasize the contrast between object and background. First, the blurred image is transformed to fit a logarithmic scale. Separately, the Scharr edge detection algorithm is applied to the blurred image. Finally, the two resulting images are merged to complete the step.

## CORE

### **0.1.5 | Global + local thresholding *(AICSSeg – MO)***

This step applies the Allen Cell & Structure Segmenter *“Masked Object Threshold” (MO)* to the output of the previous step *(0.1.4)*. First, a global threshold is performed on the image’s intensity values. Next, an Otsu threshold is applied locally for each connected component resulting from the global threshold. This creates an aggregate cellmask segmentation.

Parameters:

-	`global_method`: The type of thresholding method to be performed globally
-	`cutoff_size`: The minimum object size to advance to the local thresholding step
-	`local_adjust`: The ratio applied to the local Otsu threshold

❗ If you do not wish to adjust the local threshold set local_adjust to 1.


## POSTPROCESSING

### **0.1.6 | Remove small holes and objects**

Small holes and objects are removed to clean up the segmentation.

Parameters:

-	`hole_min`: the minimum hole width to be filled in the segmentation 
-	`hole_max`: the maximum hole width to be filled in the segmentation 
-	`min_size`: the minimum width of the objects to remain in the segmentation
-	`method`: determines if both filters are applied in 3D or on a 2D level *(slice by slice)*

❗ If the method drop-down menu is set to *“slice_by_slice”* the min_size, hole_min, and hole_max values are squared, if set to *“3D”* the values are cubed.


## POSTPOSTPROCESSING

### **0.1.7 | Select one cellmask/nuclei based on signal**

This step applies a water-shedding method to the aggregate cell mask segmentation with the nuclei objects as the sources. Each nucleus is given a unique label; the label with the largest total signal is then selected as the main cell’s cellmask. The source for said cellmask is selected as the main cell’s nucleus, and extracellular voxels are discarded.

Parameters:

-	`interior_labels_only`: whether to only consider labels within the nuclei objects
-	`watershed_method`: whether to apply the water-shedding algorithm in 3D or in 2D *(‘slice-by-slice’)*

### **0.1.8 | Segment cytoplasm**

In the penultimate step, the cytoplasm mask is created via a logical exclusive or *(XOR)* of the cellmask and the nucleus mask.

Parameter:

-	`erode_nuclei`: whether to erode the nucleus mask during the exclusion process *(makes the nucleus gap in the cytoplasm appear smaller)*


## EXPORT

### **0.1.9 | stack masks**

The final step stacks the three masks in order of: nucleus, cellmask, and cytoplasm mask. The three masks are stored in a four-dimensional format and can be viewed individually in the output image layer using the sliders in the Napari UI. Toggling 3D view is recommended for inspection purposes.


### The masks segmentation workflow is now complete. The final layer can be saved using File > Save.

### Batch processing: If you wish to utilize these settings on more than one image, or revisit the settings later, use the *“Save Workflow”* button at the bottom of the workflow editor. Be sure the settings file ends in *“masks”*.