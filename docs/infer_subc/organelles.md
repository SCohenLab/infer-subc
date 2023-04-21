# organelles sub-module

The code in this sub-module have the routines to perform a (hopefully) robust inference of subcellular objects:

+ 1️⃣-***cellmask***
+ 2️⃣-***nuclei***
+ 3️⃣-***cytoplasm***
+ 4️⃣-***lysosome***
+ 5️⃣-***mitochondria***
+ 6️⃣-***golgi***
+ 7️⃣-***peroxisome***
+ 8️⃣-***endoplasmic reticulum***
+ 9️⃣-***lipid body***

From the results of the  1️⃣-***cellmask***, 2️⃣-***nuclei***, the  3️⃣-***cytoplasm***, a mask of the cytoplasm for each cell of interest is derived.

By "inference of sub-cellular objects" we mean assigning each pixel to belonging to an organell as estimated from the florescence image in the apropriate channel.  This is done here by image processing and thresholding.

> insert an overview of our overall pipeline framework here..
> i.e. GOAL SETTING, DATA CREATION, IMAGE PROCESSING, and QUANTIFICATION
--------------
## PIPELINE OVERVIEW
### 1. GOAL SETTING ✍

#### GOAL:  Infer sub-cellular components in order to understand interactome 

To measure shape, position, size, and interaction of eight organelles/cellular components (Nuclei (NU), Lysosomes (LS),Mitochondria (MT), Golgi (GL), Peroxisomes (PO), Endoplasmic Reticulum (ER), Lipid Droplet (LD), and CELLMASK) during differentiation of iPSCs, in order to understand the Interactome / Spatiotemporal coordination.

#### summary of _OBJECTIVES_ ✅
Robust inference of subcellular objects:
- 1️⃣-***cellmask***
- 2️⃣-***nuclei***
- 3️⃣-***cytoplasm***
- 4️⃣-***lysosome***
- 5️⃣-***mitochondria***
- 6️⃣-***golgi***
- 7️⃣-***peroxisome***
- 8️⃣-***endoplasmic reticulum***
- 9️⃣-***lipid body***

### 2. DATA CREATION
> METHODS:📚📚
> 
> iPSC lines prepared and visualized on Zeiss Microscopes. 32 channel multispectral images collected.  Linear Unmixing in  ZEN Blue software with target emission spectra yields 8 channel image outputs.  Channels correspond to: Nuclei (NU), Lysosomes (LS),Mitochondria (MT), Golgi (GL), Peroxisomes (PO), Endoplasmic Reticulum (ER), Lipid Droplet (LD), and a “residual” signal.

> Meta-DATA 🏺 (artifacts)
>  - Microcope settings
>  - OME scheme
> - Experimenter observations
> - Sample, bio-replicate, image numbers, condition values, etc
>  - Dates
>  - File structure, naming conventions
>  - etc.

### 3. IMAGE PROCESSING  ⚙️🩻🔬
#### INFERENCE OF SUB-CELLULAR OBJECTS
The imported images have already been pre-processed to transform the 32 channel spectral measuremnts into "linearly unmixed" images which estimate independently labeled sub-cellular components.  Thes 7 channels (plus a residual "non-linear" signal) will be used to infer the shapes and extents of these sub-cellular components.   
A single "optimal" Z slice is chosen for each image for subsequent 2D analysis.
We will perform computational image analysis on the pictures to _segment_ (or _binarize_) the components of interest for measurement.  In other procedures we can used these labels as "ground truth" labels to train machine learning models to automatically perform the inference of these objects.
Pseudo-independent processing of the imported multi-channel image to acheive each of the 9 objecives stated above.  i.e. infering: NUCLEI, CELLMASK, CYTOPLASM, LYSOSOME, MITOCHONDRIA, GOLGI COMPLEX, PEROZISOMES, ENDOPLASMIC RETICULUM, and LIPID BODIES

### General flow for infering objects via segmentation
- (extraction) 
- Pre-processing 🌒
- Core-processing (thresholding) 🌕
- Post-processing  🌘
- (post-postprocessing) 

### QC 🚧 WIP 🚧 



### 4. QUANTIFICATION 📏📐🧮

SUBCELLULAR COMPONENT METRICS
-  extent 
-  size
-  shape
-  position

### NOTE: PIPELINE TOOL AND DESIGN CHOICES?
We want to leverage the Allen Cell & Structure Segmenter.  It has been wrapped as a [napari-plugin](https://www.napari-hub.org/plugins/napari-allencell-segmenter) but fore the workflow we are proving out here we will want to call the `aicssegmentation` [package](https://github.com/AllenCell/aics-segmentation) directly.

#### ​The Allen Cell & Structure Segmenter 
​The Allen Cell & Structure Segmenter is a Python-based open source toolkit developed at the Allen Institute for Cell Science for 3D segmentation of intracellular structures in fluorescence microscope images. This toolkit brings together classic image segmentation and iterative deep learning workflows first to generate initial high-quality 3D intracellular structure segmentations and then to easily curate these results to generate the ground truths for building robust and accurate deep learning models. The toolkit takes advantage of the high replicate 3D live cell image data collected at the Allen Institute for Cell Science of over 30 endogenous fluorescently tagged human induced pluripotent stem cell (hiPSC) lines. Each cell line represents a different intracellular structure with one or more distinct localization patterns within undifferentiated hiPS cells and hiPSC-derived cardiomyocytes.

More details about Segmenter can be found at https://allencell.org/segmenter
In order to leverage the A

## `organelles` submodule

::: infer_subc.organelles