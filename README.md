**MRI pre-preocessing for q-bold maps**



**Software Requirements**

- MATLAB (with Image Processing and Statistical Toolboxes recommended).
- SPM (Statistical Parametric Mapping) for neuroimaging analysis --> https://www.fil.ion.ucl.ac.uk/spm/.
- FSL (FMRIB Software Library) for preprocessing --> https://fsl.fmrib.ox.ac.uk/fsl/docs/#/.
Ensure FSLDIR is set correctly in the environment.
- HD-BET (Brain Extraction Tool using deep learning) --> https://github.com/MIC-DKFZ/HD-BET.
- Nifti_Util functions --> https://www.mathworks.com/matlabcentral/fileexchange/8797-tools-for-nifti-and-analyze-image

**Folders Organization**
|
├── anat/
│   ├── sub-folder/
│   │   ├── T1w/
│   │   │   └── T1w_orig.nii
│   │   ├── FLAIR/
│   │       └── FLAIR.nii
│
├── T2echoes/
│   └── T2map_e1.nii ... T2map_e8.nii
│
├── T2starechoes/
│   └── T2star_e01.nii ... T2star_e10.nii
│
├── perfusion/
│   ├── CBF.nii
│   ├── CBV.nii
│   └── dsc/
│       └── DSC.nii
│
├── T1segmentation/
│   └── 3DTumor_FLAIR.nii (In T1 space)
│
├── diffusion/
│   ├── ADC.nii
│   └── DWI.nii
│
└── DCE/
│   ├── VP.nii
│   └── DCE.nii
│
└── cerebrale_mac/
    └── Metionina_SUV.nii
