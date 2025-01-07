**MRI pre-preocessing for q-bold maps**

This repository contains code that is part of a project aimed at evaluating an innovative PET and MRI approach for assessing hypoxia, perfusion, and tissue diffusion in HGGs and deriving a combined map for the clustering of intra-tumor heterogeneity (funded by the Ministero della Salute - GR-2018-12365670).

This repository includes the main script _'Main_qbold_loop_new_met_HDbet'_ and functions to pre-process structural images (T1w and FLAIR), ADC DWI-derived maps, CBF and CBV DSC-derived maps, and VP DCE-derived maps.

It also elaborates the T2* and T2 maps by using a monoexponential fitting on T2* and T2 echoes.

The script calculates the q-BOLD maps: OEF, PO2, and CMRO2.

<p align="center">
<img src="https://github.com/NicoloPecco/MRI-Pre-processing-for-Q-bold-maps/blob/main/Figures/Panel.png" width="980" height="400">
</p>

**Note:**
The link with all data preprocess data and scripts for **Spatial Habitat imaging** will be uploaded asap.
You are free to take any piece of code you need for your anlysis!

**Software Requirements**

- MATLAB (with Image Processing and Statistical Toolboxes recommended).
- [SPM](https://www.fil.ion.ucl.ac.uk/spm/) (Statistical Parametric Mapping) for neuroimaging analysis.
- [FSL](https://fsl.fmrib.ox.ac.uk/fsl/docs/#/) (FMRIB Software Library) for preprocessing .
Ensure FSLDIR is set correctly in the environment.
- [HD-BET](https://github.com/MIC-DKFZ/HD-BET) (Brain Extraction Tool using deep learning).
- [Nifti_Util functions](https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image).

**Folders Organization**

<p align="left">
<img src="https://github.com/NicoloPecco/MRI-Pre-processing-for-Q-bold-maps/blob/main/Figures/Folders_organization.png" width="580" height="700">
</p>

**Citation**

Bailo, M., Pecco, N., Callea, M., Scifo, P., Gagliardi, F., Presotto, L., ... & Castellano, A. (2022). Decoding the heterogeneity of malignant gliomas by PET and MRI for spatial habitat analysis of hypoxia, perfusion, and diffusion imaging: a preliminary study. Frontiers in Neuroscience, 16, 885291.

Other to be updated...
