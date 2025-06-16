# **Sample Data Folder** 📂

We've included four sample images to use as input for **`infer-subc`**. The first two (neuron_1 and astrocyte) images and outputs are from experiments featured in this preprint [Neurons and astrocytes have distinct organelle signatures and responses to stress](https://www.biorxiv.org/content/10.1101/2024.10.30.621066v1). In this paper infer-subc was used to analyze differences between neurons and astrocytes at baseline and in response to stress.

To access and use the sample data in the notebooks or Napari,
follow the instructions below:

**1.** Clone the **infer-subc** repository

**2.** Follow the steps in [notebook 1.0_image_setup](/notebooks/part_1_segmentation_workflows/1.0_image_setup.ipynb) to download and correctly store the raw images in the [**sample_data**](../sample_data/) folder.

**3.** Use workflow 1.1 to segment the masks for each example image. This can be done in the notebooks or Napari.

- For the example **neuron (1)** use workflow [1.1a](/notebooks/part_1_segmentation_workflows/1.1a_infer_masks_from-composite_single_cell.ipynb)
- For the example **astrocyte** use workflow [1.1b](/notebooks/part_1_segmentation_workflows/1.1b_infer_masks_from-composite_multiple-cells.ipynb)
- For the example **neuron (2)** use workflow [1.1c](/notebooks/part_1_segmentation_workflows/1.1c_infer_masks_from-composite_neuron_with_pm.ipynb)
- For the example **ipsc** use workflow [1.1d](/notebooks/part_1_segmentation_workflows/1.1d_infer_masks_from-composite_ipsc.ipynb)

**4.** Then run workflows 1.2-1.7 twice (once per cell type).

- In the notebooks, you can switch between cell types by setting the `sample_data_type` variable equal to **"neuron_1"**, **"astrocyte"**, **"neuron_2"** or **"ipsc"**.

**5.** Carry out quantification of the sample neuron and astrocyte in part 2 using the notebooks.

##### You may also use the sample data via the [organelle-segmenter-plugin](https://github.com/ndcn/organelle-segmenter-plugin) in Napari. The raw images and workflow settings are all present within the contents of the sample data folder. For instructions on how to use the organelle-segmenter-plugin, look under the Option A section in the **[ReadMe](/README.md)**.



## Contents of Sample Data Folder 🗂️
###### This section details the expected contents of the sample data folder after the infer-subc notebooks and/or Napari plugin are run as described above. A link to the expected output is included when applicable. Of note, these links are just a reference to the expected outocome if needed. There is no need to download them because ,the notebooks in part 1 and part 2 will create them as intended. If you wish to download the original files, this can be done simply by clicking the links.

### example_astrocyte
###### This folder contains all of the files related to the segmentation of the example astrocyte

- #### raw
    ###### The intended location of the downloaded astrocyte from bioimage archive
    - [05052022_astro_control_2_Linear unmixing_0_cmle.ome.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_deconvolution/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome.tiff) (💾 1.14 GB)
- #### seg
    ###### This folder contains the resulting segmentations created in notebook 1.1b as well as the segmentation notebooks 1.2-1.7.
    - [05052022_astro_control_2_Linear unmixing_0_cmle.ome-ER.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-ER.tiff) (💾 102.6 MB)
    - [05052022_astro_control_2_Linear unmixing_0_cmle.ome-golgi.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-golgi.tiff) (💾 102.6 MB)
    - [05052022_astro_control_2_Linear unmixing_0_cmle.ome-LD.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-LD.tiff) (💾 102.6 MB)
    - [05052022_astro_control_2_Linear unmixing_0_cmle.ome-lyso.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-lyso.tiff) (💾 102.6 MB)
    - 05052022_astro_control_2_Linear unmixing_0_cmle.ome-masks_B.tiff (originally seperated into nucleus and cellmask)
        - [nucleus](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-nuc.tiff) (💾 205.2 MB)
        - [cellmask](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-cell.tiff) (💾 205.2 MB)
    - [05052022_astro_control_2_Linear unmixing_0_cmle.ome-mito.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-mito.tiff) (💾 102.6 MB)
    - [05052022_astro_control_2_Linear unmixing_0_cmle.ome-perox.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20astrocyte%20data/Z-stacks/replicate2_VD-0505/VD-0505_segmentations/05052022_astro_control_2_Linear%20unmixing_0_cmle.ome-perox.tiff) (💾 102.6 MB)
    

- #### settings
    ###### The settings (in JSON file type) used to create each of the corresponding segmentations in the seg folder. These settings can be used in napari to automatically apply the workflow settings. Each of the following files are already included in the repository, thus there is no need to download them.
    
    - ER.json
    - golgi.json
    - LD.json
    - lyso.json
    - masks_B.json
    - mito.json
    - perox.json

### example_neuron_1
###### This folder contains all of the files related to the segmentation of the example neuron

- #### raw
    ###### The intended location of the downloaded neuron from bioimage archive
    - [20230727_C2-121_unconditioned_well 10_cell 4_25nM TG_Linear unmixing_0_cmle.ome.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_deconvolution/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome.tiff) (💾 2.16 GB)
- #### seg
    ###### This folder contains the resulting segmentations created in notebook 1.1a as well as the segmentation notebooks 1.2-1.7.
    - [20230727_C2-121_conditioned_well 4_cell 3_untreated_Linear unmixing_0_cmle.ome-ER.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-ER.tiff) (💾 193.8 MB)
    - [20230727_C2-121_conditioned_well 4_cell 3_untreated_Linear unmixing_0_cmle.ome-golgi.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-golgi.tiff) (💾 193.8 MB)
    - [20230727_C2-121_conditioned_well 4_cell 3_untreated_Linear unmixing_0_cmle.ome-LD.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-LD.tiff) (💾 193.8 MB)
    - [20230727_C2-121_conditioned_well 4_cell 3_untreated_Linear unmixing_0_cmle.ome-lyso.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-lyso.tiff) (💾 193.8 MB)
    - 20230727_C2-121_conditioned_well 4_cell 3_untreated_Linear unmixing_0_cmle.ome-masks_A.tiff (originally seperated into nucleus and cellmask)
        - [nucleus](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-nuc.tiff) (💾 193.8 MB)
        - [cellmask](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-cell.tiff) (💾 193.8 MB)
    - [20230727_C2-121_conditioned_well 4_cell 3_untreated_Linear unmixing_0_cmle.ome-mito.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-mito.tiff) (💾 193.8 MB)
    - [20230727_C2-121_conditioned_well 4_cell 3_untreated_Linear unmixing_0_cmle.ome-perox.tiff](https://www.ebi.ac.uk/biostudies/files/S-BIAD1445/Neuron%20and%20astrocyte%20organelle%20signatures%20dataset/Primary%20rat%20neuron%20data/Z-stacks/replicate3_C2-121/C2-121_segmentations/20230727_C2-121_conditioned_well%204_cell%203_untreated_Linear%20unmixing_0_cmle.ome-perox.tiff) (💾 193.8 MB)

- #### settings
    ###### The settings (in JSON file type) used to create each of the corresponding segmentations in the seg folder. These settings can be used in napari to automatically apply the workflow settings. Each of the following files are already included in the repository, thus there is no need to download them.

    - ER.json
    - golgi.json
    - LD.json
    - lyso.json
    - masks_A.json
    - mito.json
    - perox.json

TBA (example_neuron_2 and example_ipsc links)
