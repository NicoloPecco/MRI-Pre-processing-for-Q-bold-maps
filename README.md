**MRI pre-preocessing for q-bold maps**

This repository contains code that are part of a project which aim to evaluate an innovative PET and MRI approach for assessing hypoxia, perfusion, and tissue diffusion in HGGs and derive a combined map for clustering of intra-tumor heterogeneity (Founded by Ministero della Salute - GR-2018-12365670). 

This repository includes the main script **'Main_qbold_loop_new_met_HDbet'** and functions to pre-process structural images (T1w and FLAIR), ADC DWI-derived maps, CBF and CBV DSC-derived maps and VP DCE-derived maps.

It also elaborate the T2* and T2 maps by using a monoexponential fitting on T2* and T2 echoes.

The script calculates the q-BOLD maps: OEF, PO2 and CMRO2.

**Note**
The link with all data and scripts for Spatial Habitat imaging will be uploaded asap.
Yuo are free to take any piece of code you need for your anlysis!

**Software Requirements**

- MATLAB (with Image Processing and Statistical Toolboxes recommended).
- SPM (Statistical Parametric Mapping) for neuroimaging analysis --> https://www.fil.ion.ucl.ac.uk/spm/.
- FSL (FMRIB Software Library) for preprocessing --> https://fsl.fmrib.ox.ac.uk/fsl/docs/#/.
Ensure FSLDIR is set correctly in the environment.
- HD-BET (Brain Extraction Tool using deep learning) --> https://github.com/MIC-DKFZ/HD-BET.
- Nifti_Util functions --> https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

**Folders Organization**

<p align="left">
<img src="https://github.com/NicoloPecco/MRI-Pre-processing-for-Q-bold-maps/blob/main/Figures/Folders_organization.png" width="500" height="700">
</p>

**Citation**

To be uploaded
