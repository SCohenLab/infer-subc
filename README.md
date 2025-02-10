README
# infer-subc
### A Python-based image analysis tool to segment and quantify the morphology, interactions, and distribution of organelles.

<img src="infer_subc\assets\README.png" width="800">
<p>

# 📒 About this project

### `infer-subc` 
- aims to create a <ins>simple</ins>, <ins>extensible</ins>, and <ins>reproducible</ins> pipeline to segment (or infer) and quantify the size, shape, interaction, and subcellular distribution of multiple intracellular organelles from confocal microscopy 🔬 images. 
- is <ins>modular</ins> 🔢 to support a variety of organelle-focused research questions. 
- can be <ins>applied broadly</ins> to many types of *in vitro* 🧫 and *in vivo* models 🐁🧬 to better understand the spatial coordination and interactome of organelles during key biological processes or disease. 

# Getting Started
## ⚙️ Setup 
`infer-subc` and the companion segmentation plugin `organelle-segmenter-plugin` for [Napari](https://napari.org/) are available via `pip`: 

```
pip install infer-subc-main
pip install infer-subc-plugin
```

A full list of dependencies and recommended setup steps are included in [env_create.sh](./env_create.sh).

## 🖍️ Part 1 - Segmentation Workflows 

> ***NOTE**: Proceed to the Organelle Quantification section below if you have already created instance segmentations using a separate method.*

The starting point for the `infer-subc` analysis pipeline is to perform instance segmentation on single or multichannel confocal microscopy images, where each channel labels a different intracellular organelle (or structure). In the infer-subc segmentation workflows, each organelle will be segmented from a single intensity channel from the microscopy image. To carry out single-cell analysis in Part 2 – Organelle Quantification, the cell area (or "mask") should be segmented as well.

> ### Compatible Organelles and Subcellular Regions 🔓🗝️
> 
> - `Lysosomes`
> - `Mitochondria`
> - `Golgi`
> - `Peroxisomes`
> - `Endoplasmic reticulum` 
> - `Lipid droplets`
> - `Cell`
> - `Nucleus`
>
>  *Outside segmentation methods can also be used to incorporate additional organelles.*

We recommend our `infer-subc` implementation for Napari called [`organelle-segmenter-plugin`](https://github.com/ndcn/organelle-segmenter-plugin) for image segmentation. This allows users to test segmentation settings for each organelle systematically, then batch process the segmentation of all organelles of interest across multiple cells using predetermined settings. Alternatively, the included set of Jupyter Notebooks can be used to work through the segmentation process step by step using functions included in the `infer-subc` module. 

If desired, alternative segmentation methods (e.g., CellProfiler, Imaris, ImageJ, etc.) can be used as the input for `Part 2 – Organelle Quantification` below.

> ### Input image formatting:
>
> We have tested the following file formats as input in both the Napari plugin and the Jupyter notebooks:
> 
> - Single or multi-channel ".tiff"/".tif" or ".czi" files
> - 2D (single Z-plane) or 3D (Z-stack) images
> - Dimension order: CZYX

### <ins>Segmentation Option A:</ins> [Napari Plugin](https://github.com/ndcn/organelle-segmenter-plugin) 🔌

The `organelle-segmenter-plugin` package is required for this method (see setup instructions above). You must segment at least one organelle and the cell mask for organelle quantification in Part 2 below.

1. Open Napari. Then drag-and-drop or use the `File` > `Open File(s)...` controls to open a single- or multi-channel confocal microscopy image. This image will be used to test the segmentation settings you want to apply during batch processing.
2. Start the plugin by navigating to `Plugin` > `Infer sub-Cellular Object Npe2 plugin` > `Workflow editor`. The plugin settings will appear as a new right-side panel.
3. In the Workflow editor, select the image you uploaded from the dropdown list. 
4. Select the workflow corresponding to your first desired organelle or the masks.
5. Adjust the parameters for each step. Use the output of each step to adjust the settings before continuing to the next step. After proceeding to a subsequent step, you cannot return to a previous step. If you need to return to a previous step, you must restart the workflow by pressing `Close Workflow` at the bottom of the panel and begin again. Your settings will not be saved; follow the next step or note down your preferred settings before closing the workflow.
6. Save the workflow settings that are compatible with your image by using the `Save Workflow` option at the bottom of the panel. *IMPORTANT: the file name should end with the same name as the workflow you are working on.*
    > 
    > <ins>**Naming Examples**</ins>: 
    >
    > For settings saved from the 0.2.lyso workflow, the following names are **acceptable**:
    > - "20241031_lyso.json"
    > - "iPSCs_lyso.json"
    > - "lyso.json"
    > <p>
    > 
    > Do **NOT** use names like:
    > - "lysosomes.json"
    > - "LS.json"
7. Close the workflow and repeat the steps above for any additional organelles and/or the masks. Save each of the workflow setting files together in the same folder.
8. Once all the settings are saved, open the batch processor by going to `Plugins` > `Infer sub-Cellular Object Npe2 plugin` > `Batch processing`. A new right-side panel will appear.
9. Load the saved workflow settings and specify the input (confocal microscopy images) and output (desired location for segmentation files to be saved) folders.
10. Click `Run`. A progress bar will allow you to track your processing.

Continue to the Quality Check section explained below BEFORE moving on to Part 2 – Organelle Quantification. 

### <ins>Segmentation Option B:</ins> [Jupyter Notebooks](/docs/nbs/overview.md) 📚
The primary purpose of the Jupyter notebooks is to walk step-by-step through each of the segmentation workflows. We hope these notebooks provide a more easily accessible resource for those who are new to Python image analysis or a more flexible platform for customization.

*The notebooks below include steps to segment a single image. A batch-processing workflow has not yet been created.*

**Step 1️⃣: Identify a single cell of interest**

Use one of the following notebooks to segment (or infer) the `cellmask` and `nuclei` from your image. Each workflow differs based on the type of fluorescent labels used and the number of cells per field of view (FOV):

- Fluorescently labeled nuclei (no cell or plasma membrane marker)
    - [one or more cells per FOV](/notebooks/part_1_segmentation_workflows/1.1_infer_masks_from-composite_with_nuc.ipynb)
- No cell, plasma membrane, or nuclei labels with
    - [one cell per field of view FOV](/notebooks/part_1_segmentation_workflows/1.1a_infer_masks_from-composite_single_cell.ipynb)
    - [one or more cells per FOV](/notebooks/part_1_segmentation_workflows/1.1b_infer_masks_from-composite_multiple-cells.ipynb)

**Step 2️⃣: Segment organelles**

Each of the organelles you wish to include in your analysis should be segmented from a single fluorescently labeled structure. Use the following notebooks to segment each organelle from a single intensity channel in the input image.

2. Infer [`lysosomes`](/notebooks/part_1_segmentation_workflows/1.2_infer_lysosome.ipynb)
3. Infer [`mitochondria`](/notebooks/part_1_segmentation_workflows/1.3_infer_mitochondria.ipynb)
4. Infer [`golgi`](/notebooks/part_1_segmentation_workflows/1.4_infer_golgi.ipynb)
5. Infer [`peroxisomes`](/notebooks/part_1_segmentation_workflows/1.5_infer_peroxisome.ipynb)
6. Infer [`endoplasmic reticulum (ER)`](/notebooks/part_1_segmentation_workflows/1.6_infer_ER.ipynb)
7. Infer [`lipid droplets`](/notebooks/part_1_segmentation_workflows/1.7_infer_lipid_droplet.ipynb)

### <ins>Quality Check and Mask Separation:</ins> [Validate segmentation results]()🔎
After segmenting all the cells in your dataset, we recommend you quality check your segmentation results by visually inspecting the images. The Segmentation Validation pipeline is included in the [Full Quantification Pipeline Notebook](/notebooks/part_2_quantification/full_quantification_pipeline.ipynb) to streamline the validation process.

**REQUIRED:** This notebook will also separate the `masks` segmentation output into separate `cell` and `nuc` (i.e., nucleus) segmentation files. This is ***required*** for Part 2 – Organelle Quantification.

🚧 *In a future version, this notebook will also include quality checks for assumptions made during quantification (i.e., only one nucleus and ER per cell, etc.).*

> ### Segmentation output formatting: 
> 
> Segmentation outputs from the Napari plugin during batch processing or the notebooks will be saved as ".tiff" files. All organelle segmentations will include a single channel. The "masks" (e.g., cell, cytoplasm, nucleus) file from the Napari plugin will stacked into a multichannel image. It must be separated into “cell” and “nuc” files for quantification (see above).

## 🧮📐 Organelle Quantification 

After all of the organelles of interest are segmented, single or multi-organelle analysis can be carried out using Jupyter Notebook-based pipeline(s). Each of the following analysis types is modular and can be used in combination or separately.

**Combined “Organelle Signature Analysis” pipeline:** 
- [Full Quantification Pipeline](./notebooks\part_2_quantification\full_quantification_pipeline.ipynb):  This notebook carries out quantification of the `morphology`, `interactions`, and `distribution` of two or more organelles within a specified region (e.g., the cell). This pipeline incorporates batch quantification for all files from a single experiment (contained in one folder) and then summarizes the quantification outputs across multiple experimental replicates. 

**Individual analysis pipelines:**

The following notebooks primarily act as a step-by-step guide to understanding each measurement type. However, they can also be used to quantify features of single organelles or pairs of organelles (interactions) from individual cells.  
- [Organelle morphology](./notebooks\part_2_quantification\1.1_organelle_morphology.ipynb)
- [Pairwise organelle interactions](./notebooks\part_2_quantification\1.2_organelle_interactions.ipynb)
- [Subcellular distribution](./notebooks\part_2_quantification\1.3_distribution_measurements.ipynb)
- [Cell/nucleus morphology](./notebooks\part_2_quantification\1.4_cell_region_morphology.ipynb)
- COMBINED ANAYSIS ONLY: [batch processing](./notebooks\part_2_quantification\1.5_combined_and_batch_processing.ipynb)
- COMBINED ANAYSIS ONLY: [per-cell summary](./notebooks\part_2_quantification\1.6_summary_stats.ipynb)

🚧 *Future implementations of these notebooks will include batch processing capabilities (e.g., multiple cells, multiple organelles) for each quantification type separately.*

# Additional Information
## Built With
A quick note on the tools and resources used...

- [`napari-allencell-segmenter`](https://github.com/AllenCell/napari-allencell-segmenter) -- We are leveraging the framework of the `napari-allencell-segmenter` plugin, which enables powerful 3D image segmentation while taking advantage of the `napari` graphical user interface. 
- [`aicssegmentation`](https://github.com/AllenCell/aics-segmentation) -- We call the `aicssegmentation` package directly to access their advanced segmentation functions.
- [`napari`](https://napari.org/stable/) -- Used as the visualization framework, a fast, interactive, multi-domensional image viewer for Python.
- [`scipy`](https://scipy.org/install/) -- Image analysis
- [`scikit-image`](https://scikit-image.org/) -- Image analysis
- [`itk`](https://itkpythonpackage.readthedocs.io/en/master/Quick_start_guide.html) -- Image analysis
- [`numpy`](https://numpy.org/) -- Under the hood computation
- [`pandas`](https://pandas.pydata.org/) -- Quantitative data manipulation

### Segmentation workflow & Napari plugin design:
Early in the development of infer-subc, we chose to leverage methods created in the `Allen Cell & Structure Segmenter` and [`napari plugin`](https://www.napari-hub.org/plugins/napari-allencell-segmenter). Although the logic of our **multi-channel** organelle segmentations required us to fork and modify their code, we hope it provides a stable but evolving base that will help manage the accumulation of technical debt. In addition to the overall logic, we particularly leverage their *workflow* paradigm, which is integral in the use of the napari plugin interface. Implementation of `infer-subc` as a Napari plugin using this framework is called [`organelle-segmenter-plugin`](https://github.com/ndcn/organelle-segmenter-plugin).

## Issues
If you encounter any problems, please file an issue with a detailed description.

## Development
Read the [CONTRIBUTING.md](CONTRIBUTING.md) file.

## License
Distributed under the terms of the [BSD-3] license.

`infer-subc` and `organelle-segmenter-plugin` are free and open-source software.

## Support of this project includes:
- [National Institutes of Health](https://www.nih.gov/) under awards T32 NS007431, F31
842 AG079622, R01NS105981, and R35GM133460
- [CZI Neurodegeneration Challenge Network (NDCN)](https://chanzuckerberg.com/science/programs-resources/neurodegeneration-challenge/)

# Publications

`infer-subc` analysis has been featured in:
1. Shannon N. Rhoads, Weizhen Dong, Chih-Hsuan Hsu, Ngudiankama R. Mfulama, Joey V. Ragusa, Michael Ye, Andy Henrie, Maria Clara Zanellati, Graham H. Diering, Todd J. Cohen, Sarah Cohen. *Neurons and astrocytes have distinct organelle signatures and responses to stress.* bioRxiv 2024.10.30.621066; doi: https://doi.org/10.1101/2024.10.30.621066

