
# infer-subc

![GitHub License](https://img.shields.io/github/license/SCohenLab/infer-subc)
![PyPI - Downloads](https://img.shields.io/pypi/dm/infer-subc?color=purple)
### A Python-based image analysis tool to segment and quantify the morphology, interactions, and distribution of organelles.

<img src="infer_subc\assets\README.png" width="800">
<p>

# 📒 About this project

### `infer-subc` 
- aims to create a <ins>reproducible</ins> pipeline to segment (or "infer") and quantify the size, shape, interaction, and subcellular distribution of multiple intracellular organelles from confocal microscopy 🔬 images. 
- is <ins>modular</ins> 🔢 to support a variety of organelle-focused research questions. 
- can be <ins>applied broadly</ins> to many types of *in vitro* 🧫 and *in vivo* models 🐁🧬 to better understand the spatial coordination and interactome of organelles during key biological processes or disease. 

# Getting Started
## ⚙️ Setup 
`infer-subc` and the companion segmentation plugin `organelle-segmenter-plugin` for [Napari](https://napari.org/) are available via `PyPI`. Install the packages as follows:

```
pip install infer-subc-main
pip install infer-subc-plugin
```

We recommend installing and using these packages in a Python environment (e.g., conda). A list of setup steps are included in [env_create.sh](./env_create.sh).

### Cloning `infer-subc`:
Cloning `infer-subc` is necessary if you are going to do any of the following:

- Run segmentation or quantification using the provided `sample data`
- If you want to modify the underlying code for specific use cases

To clone this repository, use your terminal navigate to the location on your computer where you want the clone of repository to be downloaded. Then run:

```
git clone https://github.com/SCohenLab/infer-subc.git
```
## 📂 File format
> ### <ins>Input image format</ins>:
>
> We have used the following file formats as input in both the Napari plugin and the Jupyter notebooks during development and testing of `infer-subc`:
> 
> - Single or multi-channel ".tiff"/".tif" or ".czi" files
> - 3D (Z-stack) images
> - Ideal dimension order: CZYX

> ### <ins>Segmentation output format</ins>:
>
> `infer-subc` Part 1 - Segmentation Workflows will output segmentation files as follows:
> - Single channel ".tiff" files 
>   - Subcellular regions will initially be exported from batch processing as multi-channel files
>   - See the [`quality_check_segmentations`](/notebooks/part_1_segmentation_workflows/quality_check_segmentations.ipynb) for information on how to separate these into single channels
> - The original file name will be included as the stem of the file name and a unique suffix will be appended to the end of each segmentation to signify the organelle or subcellular structure that was segmented
> - Organelle segmentation (except the ER which is *always* considered one object) will contain instance segmentations where each identified object if given a unique ID number
> - Subcellular regions ('cell', 'nucleus', 'soma', 'neurites') and the ER will only include a single labeled object per output image. In the case of the ER or neurites, where there can commonly be several disconnected component, each identified component is given the same ID number and quantified as a single object.
>
> These segmentations will act the part of the input for quantification in `infer-subc` Part 2 - Organelle Quantification. If using an alternative segmentation approach, the above listed format should be followed.

> ### <ins>Quantification input format</ins>:
>
> `infer-subc` Part 2 - Organelle Quantification will use the following files as part of the input:
> - Multi-channel intensity images from which the segmentations were derived
> - Single-channel ".tiff" organelle and subcellular region (if applicable) segmentation files
>
> Additional information related to the desired quantification methods is also required. See the [Part 2 notebooks](/notebooks/part_2_quantification/) for more information on specifics.

> ### <ins>Quantification output format</ins>:
>
> `infer-subc` Part 2 - Organelle Quantification includes two rounds of quantitafication:
> 1. Quatitative feature extraction - intensity and segmentation images are used as the input and numerical data is output. This should be run *per experimental replicate*
> 2. Per subcellular region or image summarization - quantitative data is input and summarized per subregion (if a mask if used). Multiple experimental replicates of data can be combined in this step to summarize all data that will be statistically compared.
>
> The output files in both cases have the following format:
> - One ".csv" file per quantification method (e.g., morphology, distribution, etc.)
> - Each file will begin with the "dataset name", a unique identifier for each experimental replicate 

> ### <ins>Required file structure</ins>:
> We recommend use of the following file structure:
> 1. Data for each experimental replicate should be saved in a separate folder.
> 2. All segmentation data for a biological replicate should be saved in one folder. This folder would ideal be within the same parent folder as the raw data it was derived from. *We also **highly recomment** saving the workflow settings (or batch_process_segmentation Jupyter notebook) within this folder to ensure the segmentation methods are easily identifiable in the future.* 
>    - If modifications are necessary that result in additional versions of the segmentation outputs, include the updated segmentation in a separate folder (more details on this are included in the [`quality_check_segmentations`](/notebooks/part_1_segmentation_workflows/quality_check_segmentations.ipynb) notebook).
> 3. A separate folder should be included for quantification outputs. The quantification and summary statistics can be within the same folder, if desired. This folder would ideal be within the same parent folder as the raw data and segmentation files it was derived from. *We also **highly recommend** saving the quantification notebooks used to generate the quantitative data within this folder to ensure the segmentation methods are easily identifiable in the future.* 
>    - If edits were made to segmentations or the quantification settings that result in new quantitative output, include the new analysis as a separate folder.
>
> 
> **An example file structure:**
> - 📂 experiment_1
>     - 📂 raw_data
>         - 📜 date_condition1_cell1.czi
>         - 📜 date_condition2_cell1.czi
>         - 📜 ...
>     - 📂 segmentation_data
>         - 📜 date_condition1_cell1-cell.tif
>         - 📜 date_condition1_cell1-nuc.tiff
>         - 📜 date_condition1_cell1-lyso.tiff
>         - 📜 date_condition1_cell1-mito.tiff
>         - 📜 date_condition1_cell1-golgi.tif
>         - 📜 date_condition1_cell1-perox.tiff
>         - 📜 date_condition1_cell1-ER.tiff
>         - 📜 date_condition1_cell1-LD.tiff
>         - 📜 date_condition2_cell1-cell.tiff
>         - 📜 date_condition2_cell1-nuc.tiff
>         - 📜 date_condition2_cell1-lyso.tiff
>         - 📜 date_condition2_cell1-mito.tiff
>         - 📜 date_condition2_cell1-golgi.tiff
>         - 📜 date_condition2_cell1-perox.tiff
>         - 📜 date_condition2_cell1-ER.tiff
>         - 📜 date_condition2_cell1-LD.tiff
>         - 📓 batch_process_segmentations.ipynb
>         - 📜 ...
>     - 📂 quantification_output
>         - 📜 datasetname-organelle_morphology_metrics.csv
>         - 📜 datasetname-per_org_morphology_summarystats.csv
>         - 📓 2.1_organelle_morphology.ipynb
> - 📂 experiment_2
>     - 📂 raw_data
>     - 📂 segmentation_data
>     - 📂 quantification_output

> ### <ins>Output format</ins>:
> 

## 🖍️ Part 1 - Segmentation Workflows 

> ***NOTE**: Proceed to the Organelle Quantification section below if you have already created instance segmentations of organelles and/or subcellular regions you plan to include in quantification.*

The starting point for the `infer-subc` analysis pipeline is to perform instance segmentation on single or multichannel confocal microscopy images, where each channel labels a different intracellular organelle (or structure). In the infer-subc segmentation workflows included in Part 1, each organelle will be segmented from a *single* intensity channel from the input microscopy image. Subcellular regions of interest, including the cell mask, nucleus, soma, and neurites, can be segemented to include region-specific quantitative analysis in Part 2 (see more below).

> ### Compatible Organelles and Subcellular Regions 🔓🗝️
> 
> - `Lysosomes`
> - `Mitochondria`
> - `Golgi`
> - `Peroxisomes`
> - `Endoplasmic reticulum` 
> - `Lipid droplets`
> - `Cell`/`Nucleus`
> - `Soma`/`Neurites`
>
>  *Outside segmentation methods can also be used to incorporate additional organelles or subcellular regions, if desired.*

We recommend our `infer-subc` implementation for Napari called [`organelle-segmenter-plugin`](https://github.com/SCohenLab/organelle-segmenter-plugin) for image segmentation. This allows users to optimize segmentation settings for each organelle systematically, then batch process the segmentation of all organelles of interest across multiple cells using the optimized settings. We have included a set of [Jupyter Notebooks](/notebooks/part_1_segmentation_workflows/) that include the same segmentation workflows. These notebooks are a great source of information on each workflow step, and they can act as a starting point for those who wish to modify workflow to better suite their images. 


### <ins>Segmentation Option A:</ins> [Napari Plugin](https://github.com/ndcn/organelle-segmenter-plugin) 🔌

The `organelle-segmenter-plugin` package is required for this method (see setup instructions above). You must segment at least one organelle and the cell mask for organelle quantification in Part 2 below.

1. Open Napari. Then drag-and-drop or use the `File` > `Open File(s)...` controls to open a single- or multi-channel confocal microscopy image. This image will be used to test the segmentation settings you want to apply during batch processing.
2. Start the plugin by navigating to `Plugin` > `Infer sub-Cellular Object Npe2 plugin` > `Workflow editor`. The plugin settings will appear as a new right-side panel.
3. In the Workflow editor, select the image you uploaded from the dropdown list. 
4. Select the workflow corresponding to your first desired organelle or the masks.
5. Adjust the parameters for each step, one at a time. You can adjust the settings within a single step as many times as you would like; each time a step is run, a new output layer appears. After proceeding to a subsequent step, you cannot return to a previous step. If you need to return to a previous step, you must restart the workflow by pressing `Close Workflow` at the bottom of the panel and begin again. Your settings will not be saved automatically; follow the next step or note down your preferred settings before closing the workflow.
6. Once you are satisified with the workflow settings you've selected (*tip: we recommend typing them one a variety of images/experimental conditions to assess robustness and refine settings as needed*), save the workflow settings that are compatible with your image by using the `Save Workflow` option at the bottom of the panel. *IMPORTANT: the file name should end with the same name as the workflow you are working on.*
    > 
    > <ins>**Naming Examples**</ins>: 
    >
    > For settings saved from the 0.2.lyso workflow, the following names are **acceptable**:
    > - "20241031_lyso.json"
    > - "iPSCs_lyso.json"
    > - "lyso.json"
    > <p>
    > 
    > Do **NOT** use names like: (does not end in workflow suffix)
    > - "lysosomes.json" 
    > - "LS.json"
7. Close the workflow and repeat the steps above for any additional organelles and/or the masks. Save each of the workflow setting files together in the same folder.
8. Once all the settings are saved, open the batch processor plugin in Napari by going to `Plugins` > `Infer sub-Cellular Object Npe2 plugin` > `Batch processing`. A new right-side panel will appear.
9. Load the saved workflow settings (all of them can be processed at the same time) and specify the input (confocal microscopy images) and output (desired location for segmentation files to be saved) folders.
10. Click `Run`. A progress bar will allow you to track your processing.

Continue to the Quality Check section explained below **BEFORE** moving on to Part 2 – Organelle Quantification.

### <ins>Segmentation Option B:</ins> [Jupyter Notebooks](/docs/nbs/overview.md) 📚
We have supplied the same analysis methods available in the Napari plugin in Jupyter Notebook format. The primary purpose of the Jupyter notebooks is to walk step-by-step through each of the segmentation workflows, linking the underlying code to each step in the segmentation workflows. We hope these notebooks provide a more easily accessible resource for those who are new to Python image analysis or a more flexible platform for customization.

*The Jupyter Notebooks can be used in a similar fashion as the Napari plugin: 1) optimize segementation settings for each workflow; 2) batch process multiple segmentation workflows simultaneously on a set of images.*

1. Download the setup notebook (1.0) and the segmentation workflow notebooks (1.1-1.8) needed for your analysis. All segmentation workflow notebook can be found in this repository under `notebooks`>[`part_1_segmentation_workflows`](/infer-subc/notebooks/part_1_segmentation_workflows/).
2. Work through notebook [1.0_image_setup](/infer-subc/notebooks/part_1_segmentation_workflows/1.0_image_setup.ipynb) to ensure your images are compatible with the current infer-subc file readering and information extraction approaches. Any necessary updates needed for your images can be tested and implemented here. The steps presented in this notebook will be used to open raw files and read metadata in all other part 1 notebooks.
3. Use notebooks 1.1 through 1.8 to determine the appropriate segmentation settings for your images (*tip: we recommend typing segmentation settings one a variety of images/experimental conditions to assess robustness and refine settings as needed*). The settings implemented in these notebooks will be used as a reference when setting up batch processing in the next step.
4. After you have determined the optimal segmentation settings for each desired workflow, work through the [batch_process_segmentation](/infer-subc/notebooks/part_1_segmentation_workflows/batch_process_segmentations.ipynb) notebook to batch process a series of images (all from the same folder). 

Continue to the Quality Check section explained below **BEFORE** moving on to Part 2 – Organelle Quantification.

### <ins>Quality Check and Mask Separation:</ins> [Validate segmentation results]()🔎
After segmenting all the cells in your dataset, we recommend you quality check your segmentation results by visually inspecting the images. The [quality_check_segmentation](/infer-subc/notebooks/part_1_segmentation_workflows/quality_check_segmentations.ipynb) notebook walks you through the quality checking process we recommend. This notebook also separates the `masks` segmentation output into separate `cell` and `nuc` (i.e., nucleus) segmentation files as well as the optional `soma_neurite` segmentation output into separate `soma` and `neurite` segmentation files, if you are include them in your analysis. This is ***`REQUIRED`*** for Part 2 – Organelle Quantification. The notebook also ensures your data meet several assumptions necessary for quantification.


> ### Segmentation output formatting: 
> 
> Segmentation outputs from the Napari plugin or notebook during batch processing will be saved as ".tiff" files. All organelle segmentations will include a single channel. The "masks" (e.g., cell, nucleus) and "soma_neurites" files will be stacked into a multichannel image. They ***must*** be separated into “cell” and “nuc” (or "soma" and "neurites") files before quantification (see the Quality Check section above).

## 🧮📐 Organelle Quantification 

After all of the organelles of interest are segmented, single or multi-organelle analysis can be carried out using Jupyter Notebook-based pipeline(s). There are two main analysis approaches you can utilize:

**1. Individual analysis pipelines:**

The following notebooks primarily act as a step-by-step guide to understanding each measurement type. However, they can also be used to quantify features of single organelles or pairs of organelles (interactions) from individual cells.  
- [Organelle morphology](./notebooks/part_2_quantification/2.1_organelle_morphology.ipynb)
- [Organelle interactions](./notebooks/part_2_quantification/2.2_organelle_interactions.ipynb)
- [Subcellular distribution](./notebooks/part_2_quantification/2.3_organelle_distribution.ipynb)
- [Regions morphology](./notebooks/part_2_quantification/2.4_cell_region_morphology.ipynb)

**2. Combined “Organelle Signature Analysis” pipeline:** 
- [Full Quantification Pipeline](./notebooks/part_2_quantification/organelle_signature_analysis.ipynb):  This notebook combines the modular analyses into a single pipeline that quantifies the `morphology`, `interactions`, and `distribution` of two or more organelles within a specified region (e.g., the cell) or the whole image. This pipeline batch processes quantification for all files from a single experiment (contained in one folder) and then summarizes the quantification outputs across multiple experimental replicates. 

### <ins>Quantification via Jupyter Notebooks:</ins>
1. Download the setup notebook (2.0) and the quantification notebook(s) (2.1-2.4 or organelle_signature_analysis) you wish to use for your quantitative analysis. All quantification notebook can be found in this repository under `notebooks`>[`part_2_quantification`](./notebooks/part_2_quantification/).
2. Work through notebook [2.0_quantification_setup](./notebooks/part_2_quantification/2.0_quantification_setup.ipynb) to ensure your data are compatible with the current infer-subc file readering and information extraction approaches. If you utilized the segmentation workflows available in Part 1 of `infer-subc`, your setup should be straightfoward. However, any necessary updates needed for your images can be tested and implemented here. The steps presented in this notebook will be used to open raw and segmentation files in all other part 2 notebooks.
3. Use notebooks 2.1 through 2.4 or the organelle_signature_analysis notebook to quantify data from each experimental replicate, then summarize the data per region or image across multiple replicates. *See the file organization in the Part 1 section above for reference on how files should be organized.*


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

