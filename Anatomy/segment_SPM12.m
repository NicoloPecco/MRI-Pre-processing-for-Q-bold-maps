%% New Segmentation with SPM12
% Segment image resulting in up to 6 segments.
% - Stepahn Kaczmarz
% - 01JUN15

% - CHPSK10JUN15
%   + New additional Brainmask including CSF (GM & WM & CSF)
% - AK20SEP15
%   + rename 'cBrMsk...' to select brain mask with other segments


%% Main program

function segment_SPM12(path, segments_nb)
% ------------------------------------------------------------------------
% INITIALIZE
% ------------------------------------------------------------------------

path_SPM = spm('Dir'); 

% Default: All 6 segments
if (nargin == 1)
  segments_nb = 6;
  fprintf('\nChoose default value of 6 segments as no segment number was defined.\n');
end

% Check if segments number is a number
if ~isa(segments_nb, 'numeric')
  fprintf('\nInvalid selection of segment number\n');
  return;
end

% Number of Gaussians per segment
ngaus = [2, 2, 2, 3, 4 ,2];


% ------------------------------------------------------------------------
% SET UP PARAMETERS
% ------------------------------------------------------------------------

% Set path and defaults
matlabbatch{1}.spm.spatial.preproc.channel.vols     = {path};
matlabbatch{1}.spm.spatial.preproc.channel.biasreg  = 0.001;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
matlabbatch{1}.spm.spatial.preproc.channel.write    = [0 1];

% Iterate through all segments
for i=1:segments_nb
  matlabbatch{1}.spm.spatial.preproc.tissue(i).tpm      = ...
                        {fullfile(path_SPM, ['tpm/TPM.nii,' int2str(i)])};
  matlabbatch{1}.spm.spatial.preproc.tissue(i).ngaus    = ngaus(i);
  matlabbatch{1}.spm.spatial.preproc.tissue(i).native   = [1 0];
  matlabbatch{1}.spm.spatial.preproc.tissue(i).warped   = [0 0];
end

% Warp image paramters
matlabbatch{1}.spm.spatial.preproc.warp.mrf         = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup     = 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg         = [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg      = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm        = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp        = 3;
matlabbatch{1}.spm.spatial.preproc.warp.write       = [0 0];


% ------------------------------------------------------------------------
% SEGMENT
% ------------------------------------------------------------------------

% Start segmentation job
spm('defaults', 'FMRI');
spm_jobman('initcfg');
spm_jobman('serial', matlabbatch);

% ------------------------------------------------------------------------
% Brainmask
% ------------------------------------------------------------------------

% Call brainmask function
brainmask(path)

%% Generate Brainmask
% Using segmented GM, WM & CSF information

function brainmask(path)
% ------------------------------------------------------------------------
% LOAD GM, WM & CSF IMAGES
% ------------------------------------------------------------------------

% Get MPR path
[dir, ~, ~] = fileparts(path);

% Load segmented nii filenames and files
path_load_segmentated = spm_select('FPList', dir, '^c.*\.nii$');
V           = spm_vol(path_load_segmentated);

% GM, WM & CSF images
g = spm_read_vols(V(1));
w = spm_read_vols(V(2));
c = spm_read_vols(V(3));


% ------------------------------------------------------------------------
% GENERATE BRAINMASK (WM & GM)
% ------------------------------------------------------------------------

% Define brain mask
b    = w;

% Build a 3x3x3 seperable smoothing kernel
kx=[0.75 1 0.75];
ky=[0.75 1 0.75];
kz=[0.75 1 0.75];
sm=sum(kron(kron(kz,ky),kx))^(1/3);
kx=kx/sm; ky=ky/sm; kz=kz/sm;

% Erosions and conditional dilations
niter = 32;
for j=1:niter
    if j>2, th=0.15; else th=0.6; end % Dilate after two its of erosion.
    for i=1:size(b,3)
        gp = double(g(:,:,i));
        wp = double(w(:,:,i));
        bp = double(b(:,:,i));
        bp = ((bp>th).*(wp+gp));
        b(:,:,i) = uint8(round(bp));
    end
    spm_conv_vol(b,b,kx,ky,kz,-[1 1 1]);
end

th = 0.05;
for i=1:size(b,3)
    gp       = g(:,:,i);
    wp       = w(:,:,i);
    cp       = c(:,:,i);
    bp       = b(:,:,i);
    bp       = ((bp>th).*(wp+gp))>th;
    g(:,:,i) = uint8(round(gp.*bp./(gp+wp+cp+eps)));
    w(:,:,i) = uint8(round(wp.*bp./(gp+wp+cp+eps)));
    c(:,:,i) = uint8(round(cp.*bp./(gp+wp+cp+eps)+cp.*(1-bp)));
    b(:,:,i) = uint8(bp);
end


% ------------------------------------------------------------------------
% SAVE BRAINMASK (WM & GM)
% ------------------------------------------------------------------------

V_b = V(1);
V_b.fname    = fullfile( dir, 'cBrMsk.nii');
spm_write_vol( V_b, b);

% ------------------------------------------------------------------------
% EXTENDED CSF BRAINMASK (WM & GM & CSF)
% ------------------------------------------------------------------------

% Generate Brain mask incl. CSF
CSFmask = c > 0.75;
BRmsk_CSF = b | CSFmask;

% Save extended Brain Mask
V_b2 = V(1);
V_b2.fname    = fullfile( dir, 'cBrMsk_CSF.nii');
spm_write_vol( V_b2, BRmsk_CSF);
