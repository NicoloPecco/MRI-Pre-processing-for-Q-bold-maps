%% Qbold evaluation maps
clear all;close all;clc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Needed folder and file:
%       anat: 
%           subfolder:
%                T1w: T1w_orig.nii
%                FLAIR: FLAIR.nii
%
%       T2echoes: named T2map_e1.nii .. T2map_e8.nii
%
%       T2starechoes: named T2star_e01.nii .. T2star_e10.nii
%
%       Perfusion: CBF.nii, CBV.nii
%           subfolder:
%               dsc: DSC.nii
%
%       T1segmentation: 3DTumor_FLAIR.nii (In T1 space)
%       
%       Diffusion: ADC.nii DWI.nii
%       
%       DCE: VP.nii DCE.nii


setenv('FSLDIR','/usr/local/fsl');  % this to tell where FSL folder is
setenv('FSLOUTPUTTYPE', 'NIFTI_GZ'); % this to tell what the output type would be

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%     ^^^^^^^      |^^^^|  /^^^^\ |       |^^^^\     %%%%%%%%%%%%')
disp('%%%%%%%%%%%%    /       \     |    | /      \|       |     \    %%%%%%%%%%%%')
disp('%%%%%%%%%%%%   |         |    |____/|       ||       |      |   %%%%%%%%%%%%')
disp('%%%%%%%%%%%%   |         | ---|    \|       ||       |      |   %%%%%%%%%%%%')
disp('%%%%%%%%%%%%   |       \\|    |    |\       /|       |     /    %%%%%%%%%%%%')
disp('%%%%%%%%%%%%    \_______\\    |____/ \_____/ |______ |____/     %%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%% Q-Bold Evaluation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')

%% Load/Insert input:

% Correct motion in T2 and T2* Echoes
motT2(1)=1; % 1 if you want to motion correct T2 echoes
motT2(2)=1; % 1 if you want to motion correct T2 star echoes
DWI_present=1;
%Check the header file to insert the following value in ms! 
T2map_TR=1000;
T2star_TR=100;
DSC_TR=1800;
display=0;

fprintf('TEs for T2 echoes, T2*echoes and DSC are taken from nifti header.\n')
fprintf('Check if TR for T2 echoes, T2*echoes and DSC are corrected:\n')
fprintf(['T2 echoes TR: ',num2str(T2map_TR),'\n'])
fprintf(['T2*echoes: ',num2str(T2star_TR),'\n'])
fprintf(['DSC: ',num2str(DSC_TR),'\n'])

ANS1=spm_input('TRs value correct?:','-1','b','Yes|No',[1,0],1);

if ANS1
    
   fprintf('Choose the subject path with all folders and files.\n')
 
else 
    error('Uncorrect TRs. Set them correctly..')
end

%% Main Paths

fprintf('\n Creating main Path..\n')

path=uigetdir();
path=[path,'/'];

all_patients=1;

if all_patients==1

Pazienti=dir(path);
Pazienti(ismember({Pazienti.name},{'.','..','.DS_Store'}))=[];

else

Pazienti.name='SinglePazientAnalysis';

end

for PZ = 1:size(Pazienti,1)
    
    path_general=[path,Pazienti(PZ).name];
    path_general=[path_general,'/'];
    path_T1w=[path_general,'/anat/T1w'];
    path_FLAIR=[path_general,'/anat/FLAIR'];
    pathT2starEchoes=([path_general,'T2star_echoes']);
    pathT2mapEchoes=([path_general,'T2_echoes']);
    pathDSC=([path_general,'/Perfusion/dsc/']);
    pathDSCMaps=([path_general,'/Perfusion/']);
    pathdiffusion=([path_general,'/Diffusion/']);
    patht1segment=([path_general,'/T1segmentaion/']);
    pathDSC_matlab=([path_general,'/Perfusion_for_dsc']);
    pathDCE=([path_general,'DCE']);
    pathMet=([path_general,'/cerebrale_mac/']);

fprintf(['\n Working on patients: ',Pazienti(PZ).name,'...\n'])

%% Anat processing:

if isfile([path_FLAIR,'/rFLAIR_bet.nii']) 

    fprintf('\n Anatomy already Processed..\n')

else

fprintf('\n Orienting to standard FLAIR and T1...\n')

PG = spm_select('FPList', path_T1w, '^T.*.nii$');
V1 = spm_vol(PG);
dim_T1mask = V1.dim;

% Resizing T1 isotropic
fprintf('\n Resizing T1 to isotropic voxel size 1x1x1...\n')
Risizeimage(path_T1w,1)

% FSL orient to standard

command_Flair_or2std=['fslreorient2std ', path_FLAIR,'/FLAIR.nii ', path_FLAIR,'/FLAIR_reor.nii'];
system(command_Flair_or2std)
command_T1_or2std=['fslreorient2std ', path_T1w,'/T1w_orig.nii ', path_T1w,'/T1w_orig_reor.nii'];
system(command_T1_or2std)
delete([path_T1w,'/T1w_orig.nii'])
delete([path_FLAIR,'/FLAIR.nii'])
system(['gunzip ',path_FLAIR,'/*'])
system(['gunzip ',path_T1w,'/*'])
movefile([path_T1w,'/T1w_orig_reor.nii'],[path_T1w,'/T1w_orig.nii'])
movefile([path_FLAIR,'/FLAIR_reor.nii'],[path_FLAIR,'/FLAIR.nii'])

command_T1_or2std=['fslreorient2std ', patht1segment '/3DTumor_FLAIR.nii ', patht1segment '/3DTumor_FLAIR_reor.nii '];
system(command_T1_or2std)
command2=['gunzip ' patht1segment '/*'];
system(command2)
movefile([patht1segment,'/3DTumor_FLAIR.nii'],[patht1segment,'/orig_3DTumor_FLAIR.nii'])
movefile([patht1segment,'/3DTumor_FLAIR_reor.nii'],[patht1segment,'/3DTumor_FLAIR.nii'])

% Hd-bet
fprintf('\n Betting FLAIR and T1...\n')
command_Flair_bet=['hd-bet -i ', path_FLAIR,'/FLAIR.nii -device cpu -mode fast -tta 0'];
system(command_Flair_bet)
command_T1_bet=['hd-bet -i ', path_T1w,'/T1w_orig.nii -device cpu -mode fast -tta 0'];
system(command_T1_bet)
delete([path_T1w,'/T1w_orig.nii'])
delete([path_FLAIR,'/FLAIR.nii'])
delete([path_FLAIR,'/FLAIR_bet_mask.nii.gz'])
system(['gunzip ',path_FLAIR,'/*'])
system(['gunzip ',path_T1w,'/*'])
movefile([path_T1w,'/T1w_orig_bet.nii'],[path_T1w,'/T1w_orig.nii'])
movefile([path_T1w,'/T1w_orig_bet_mask.nii'],[path_T1w,'/T1_Mask.nii'])
movefile([path_FLAIR,'/FLAIR_bet.nii'],[path_FLAIR,'/FLAIR.nii'])

% Segmenting T1
fprintf('\n Segmenting T1 ...\n')
process_anatomy(path_T1w,path_FLAIR)

movefile([path_T1w,'/T1w_orig.nii'],[path_T1w,'/T1w.nii'])

% Coregister FLAIR
system(['/usr/local/fsl/bin/flirt -in ',path_FLAIR,'/FLAIR.nii -ref ',path_T1w,'/T1w.nii -out ',path_FLAIR,'/rFLAIR.nii -omat ',path_FLAIR,'/Flair_to_T1.mat']);
system(['gunzip ',path_FLAIR,'/*'])

% Checking mask:

PO1 = spm_select( 'FPList', path_FLAIR, '^F.*.nii$');
V1 = spm_vol(PO1);
dim_Flair_mask = V1.dim;

PO1 = spm_select( 'FPList', patht1segment, '^3.*FLAIR.nii$');
Vm = spm_vol(PO1);
dim_Flair_mask = Vm.dim;

if sum(dim_Flair_mask==dim_T1mask)==3

    system(['flirt -in ' patht1segment '/3DTumor_FLAIR.nii  -ref ' path_T1w, '/T1w.nii -out ',patht1segment,'/r3DTumor_FLAIR -applyisoxfm 1 -interp nearestneighbour']);
    command2=['gunzip ' patht1segment '/*'];
    system(command2)
    movefile([patht1segment,'/3DTumor_FLAIR.nii'],[patht1segment,'/previousT1_3DTumor_FLAIR.nii'])
    movefile([patht1segment,'/r3DTumor_FLAIR.nii'],[patht1segment,'/3DTumor_FLAIR.nii'])

else
    system(['/usr/local/fsl/bin/flirt -in ',patht1segment,'/3DTumor_FLAIR.nii -ref ',path_T1w,'/T1w.nii -applyxfm -init ',path_FLAIR,'/Flair_to_T1.mat',' -out ',patht1segment,'/r3DTumor_FLAIR.nii -interp nearestneighbour'])
    command2=['gunzip ' patht1segment '/*'];
    system(command2)
    movefile([patht1segment,'/3DTumor_FLAIR.nii'],[patht1segment,'/previousFLAIR_3DTumor_FLAIR.nii'])
    movefile([patht1segment,'/r3DTumor_FLAIR.nii'],[patht1segment,'/3DTumor_FLAIR.nii'])
end

% Betting T1 and FLAIR with T1 mask

movefile([path_T1w,'/T1w.nii'],[path_T1w,'/T1w_bet.nii'])
movefile([path_FLAIR,'/rFLAIR.nii'],[path_FLAIR,'/rFLAIR_bet.nii'])

WhiteStrip_Norm(path_T1w,path_FLAIR,patht1segment)

%% CHECK for FLAIR T1 and mask:

strutturaNifti(1).file = load_untouch_nii ([path_T1w,'/T1w_bet.nii']);
strutturaNifti(2).file = load_untouch_nii ([path_FLAIR,'/rFLAIR_bet.nii']);
strutturaNifti(3).file = load_untouch_nii ([patht1segment,'3DTumor_FLAIR.nii']);

disp(['T1 size: ',num2str(size(strutturaNifti(1).file.img))]);
disp(['FLAIR size: ',num2str(size(strutturaNifti(2).file.img))]);
disp(['SEG size: ',num2str(size(strutturaNifti(3).file.img))]);

multttt=(double(strutturaNifti(1).file.img)).*(double(strutturaNifti(2).file.img)).*(logical(strutturaNifti(3).file.img));

end

%% Fit T2s:

if isfile([pathT2mapEchoes,'/T2map_fin.nii']) && isfile([pathT2starEchoes,'/rT2star_uncorr_fin.nii']) 

    fprintf('\n T2 map and T2star map already processed..\n')

else

fprintf('\n Motion Correction of T2s...\n') 

[TET2,TET2s,flagT2,flagT2s]=Motion_Coregistr_T2s(pathT2mapEchoes,pathT2starEchoes,motT2);

fprintf('\n Monoexp Fitting of T2 and T2* Echoes...\n') 

THR=150;

FitMonoexpT2s(T2map_TR,T2star_TR,TET2,TET2s,pathT2mapEchoes,pathT2starEchoes,THR);

% Coreg T2star a T2 map

fprintf('\n Coregistering T2s map in T2 mapping space (Global Ref.)...\n')

Coreg_T2star1stEcho_to_T2map1stEcho(pathT2mapEchoes,pathT2starEchoes)

if display==1
spm_check_registration([pathT2mapEchoes, '/T2map_fin.nii'],[pathT2starEchoes, '/rT2star_uncorr_fin.nii']);
spm_input('Check Coregistration:','-1','b','Continue|--',[1,0]);
end

end

%% Extracting B0 from DSC for coregistration to first echo t2:

if isfile([pathDSCMaps,'/rCBV_masked.nii']) && isfile([pathDSCMaps,'rCBF_masked.nii']) 

        fprintf('\n CBV map and CBF map already processed..\n')

else

fprintf('\n Extractin B0 from DSC...\n')
   
    % Extract B0 from DSC

system(['/usr/local/fsl/bin/fslsplit ',pathDSC,'DSC.nii ',pathDSC,' -t']);
movefile([pathDSC,'/0000.nii.gz'],[pathDSC,'/B0.nii.gz'])
delete([pathDSC,'0*'])

% command_B0_DSC_bet=['hd-bet -i ', pathDSC,'/B0_reor.nii.gz -device cpu -mode fast -tta 0'];
% system(command_B0_DSC_bet)

% movefile([pathDCE,'/0000.nii.gz'],[pathDCE,'/B0.nii.gz'])

system(['/usr/local/fsl/bin/bet2 ',pathDSC,'B0.nii.gz ',pathDSC, 'B0_bet.nii.gz -m -f 0.2']);

    % Coreg CBV to T1

%fprintf('/n Correcting CFV and CBV values with ASL...')

%Correct_ASL(pathDSCMaps) % cooreggo i valori di CBF and CBV DSC

fprintf('\n Coregisterning all Perfusion to T1...\n')

% fprintf('\n Coregisterning CBV to B0...\n')
system(['/usr/local/fsl/bin/flirt -in ',pathDSCMaps,'/CBV.nii -ref ',pathDSC,'/B0.nii.gz -out ',pathDSCMaps,'/rCBV.nii -omat ',pathDSCMaps,'/CBV_to_B0.mat']);
system(['/usr/local/fsl/bin/flirt -in ',pathDSCMaps,'/CBF.nii -ref ',pathDSC,'/B0.nii.gz -out ',pathDSCMaps,'/rCBF.nii -omat ',pathDSCMaps,'/CBF_to_B0.mat']);

system(['gunzip  -f ',pathDSC,'/*'])
system(['gunzip  -f ',pathDSCMaps,'/*'])

% Coreg_CBV_to_B0(pathDSCMaps,pathDSC)

% fprintf('\n Masking CBV and CBF with B0 mask...\n')
%

command=['/usr/local/fsl/bin/fslmaths ' pathDSCMaps 'rCBV.nii -mul ' pathDSC 'B0_bet.nii.gz_mask.nii ' pathDSCMaps 'rCBV_masked.nii'];
command1=['/usr/local/fsl/bin/fslmaths ' pathDSCMaps 'rCBF.nii -mul ' pathDSC 'B0_bet.nii.gz_mask.nii ' pathDSCMaps 'rCBF_masked.nii'];

system(command);
system(command1);
system(['gunzip  -f ',pathDSCMaps,'/*'])

fprintf('\n Coregisterning CBV and B0 to T1...\n')

% Coreg_CBV_to_T1(pathDSCMaps,path_T1w)

command_B0_or2std=['fslreorient2std ', pathDSC,'/B0_bet.nii ', pathDSC,'/B0_bet_reor.nii'];
system(command_B0_or2std)

% system(['/usr/local/fsl/bin/flirt -in ',pathDSC,'/B0_bet_reor.nii.gz -ref ',path_T1w,'/T1w_bet_norm.nii -out ',pathDSC,'/rB0.nii -omat ',pathDSC,'/B0_to_T1.mat -cost normmi -searchrx -180 180 -searchry -180 180  -searchrz -180 180 -dof 12 -interp trilinear']);
system(['/usr/local/fsl/bin/flirt -in ',pathDSC,'/B0_bet_reor.nii.gz -ref ',path_T1w,'/T1w_bet_norm.nii -out ',pathDSC,'/rB0.nii -omat ',pathDSC,'/B0_to_T1.mat -dof 6 -cost mutualinfo -interp trilinear']);
system(['/usr/local/fsl/bin/flirt -in ',pathDSCMaps,'/rCBV.nii -ref ',path_T1w,'/T1w_bet_norm.nii -applyxfm -init ',pathDSC,'/B0_to_T1.mat',' -out ',pathDSCMaps,'/rrCBV_masked.nii'])
system(['/usr/local/fsl/bin/flirt -in ',pathDSCMaps,'/rCBF.nii -ref ',path_T1w,'/T1w_bet_norm.nii -applyxfm -init ',pathDSC,'/B0_to_T1.mat',' -out ',pathDSCMaps,'/rrCBF_masked.nii'])

command=['/usr/local/fsl/bin/fslmaths ' pathDSCMaps 'rrCBV_masked.nii -mul ' path_T1w '/T1_Mask.nii ' pathDSCMaps 'rCBV_masked.nii'];
command1=['/usr/local/fsl/bin/fslmaths ' pathDSCMaps 'rrCBF_masked.nii -mul ' path_T1w '/T1_Mask.nii ' pathDSCMaps 'rCBF_masked.nii'];

system(command);
system(command1);

delete([pathDSCMaps 'rCBV_masked.nii'])
delete([pathDSCMaps 'rCBF_masked.nii'])
delete([pathDSCMaps 'rrCBV_masked.nii.gz'])
delete([pathDSCMaps 'rrCBF_masked.nii.gz'])

system(['gunzip - f ',pathDSCMaps,'/*'])

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfile([pathDCE,'/rrVP_masked.nii']) && isfile([pathDCE,'/rB0_bet.nii']) 

        fprintf('\n CBV map and CBF map already processed..\n')

else
    % Extract B0 from DCE

fprintf('\n Extractin B0 from DCE...\n')

system(['/usr/local/fsl/bin/fslsplit ',pathDCE,'/DCE.nii ',pathDCE,'/',' -t']);
movefile([pathDCE,'/0000.nii.gz'],[pathDCE,'/B0.nii.gz'])
delete([pathDCE,'/0*'])
system(['/usr/local/fsl/bin/bet2 ',pathDCE,'/B0.nii.gz ',pathDCE, '/B0_bet.nii.gz -m -f 0.2']);
system(['gunzip -f ',pathDCE,'/B0_bet.nii.gz']);
fprintf('\n Coregisterning VP to B0 T1...\n')

Coreg_VP_to_B0(pathDCE) 

fprintf('\n Masking VP with B0 mask...\n')

command=['/usr/local/fsl/bin/fslmaths ' pathDCE '/rVP.nii -mul ' pathDCE '/B0_bet.nii.gz_mask.nii.gz ' pathDCE '/rVP_masked.nii'];

system(command);
system(['gunzip  -f ',pathDCE,'/*'])

fprintf('\n Coregisterning VP and B0 to T1...\n')

Coreg_VP_to_T1(pathDCE,path_T1w)

command=['/usr/local/fsl/bin/fslmaths ' pathDCE '/rrVP_masked.nii -mul ' path_T1w '/T1_Mask.nii ' pathDCE '/rrVP_masked.nii'];

system(command);
system(['gunzip  -f ',pathDCE,'/*'])

V = spm_vol([pathDCE '/rrVP_masked.nii']);
[rrVP_masked, ~] = spm_read_vols(V);
rrVP_masked(isnan(rrVP_masked))=0;
V.fname=[pathDCE,'/rrVP_masked.nii'];
spm_write_vol(V,rrVP_masked);

if display==1
spm_check_registration([path_T1w,'/T1w_bet.nii'],[pathDSCMaps,'rrCBF_masked.nii'],[pathDSCMaps,'rrCBV_masked.nii'],[pathDCE,'/rrVP_masked.nii']);
spm_input('Check Coregistration:','-1','b','Continue|--',[1,0]);
end

end

%% Extracting B0 from DWI for coregistration to first echo t2:

if isfile([pathdiffusion,'/rrADC_masked.nii']) && isfile([pathdiffusion,'rB0_bet.nii']) 

        fprintf('\n ADC map already processed..\n')

else

    if DWI_present==1
% Extract B0
fprintf('\n Extractin B0 from diffusion...\n')
system(['/usr/local/fsl/bin/fslsplit ',pathdiffusion,'DWI.nii ',pathdiffusion,' -t']);
movefile([pathdiffusion,'/0002.nii.gz'],[pathdiffusion,'/B0.nii.gz'])
delete([pathdiffusion,'0*'])
system(['/usr/local/fsl/bin/bet2 ',pathdiffusion,'B0.nii.gz ',pathdiffusion, 'B0_bet.nii.gz -m -f 0.2']);
system(['gunzip -f ',pathdiffusion,'B0_bet.nii.gz']);

% Coreg CBV to T1

fprintf('\n Coregisterning ADC to B0 T1...\n')

Coreg_ADC_to_B0(pathdiffusion)

fprintf('\n Masking ADC with B0 mask...\n')

command=['/usr/local/fsl/bin/fslmaths ' pathdiffusion '/rADC.nii -mul ' pathdiffusion '/B0_bet.nii.gz_mask.nii.gz ' pathdiffusion '/rADC_masked.nii'];

system(command);
system(['gunzip  -f ',pathdiffusion,'/*'])

fprintf('\n Coregistering ADC and B0 map to T1...\n')

Coreg_ADC_to_T1(pathdiffusion,path_T1w)

command=['/usr/local/fsl/bin/fslmaths ' pathdiffusion '/rrADC_masked.nii -mul ' path_T1w '/T1_Mask.nii ' pathdiffusion '/rrADC_masked.nii'];

system(command);
system(['gunzip  -f ',pathdiffusion,'/*'])

V = spm_vol([pathdiffusion '/rrADC_masked.nii']);
[rrADC_masked, ~] = spm_read_vols(V);
rrADC_masked(isnan(rrADC_masked))=0;
V.fname=[pathdiffusion,'/rrADC_masked.nii'];
spm_write_vol(V,rrADC_masked);


if display==1
spm_check_registration([pathdiffusion, '/rrADC_masked.nii'],[path_T1w, '/T1w_bet.nii']);
spm_input('Check Coregistration:','-1','b','Continue|--',[1,0]);
end
else
    disp('No DWI data present.')
end
end

%% Q-BOLD EVAL

if isfile([pathT2mapEchoes,'/ST2map_fin.nii']) && isfile([pathT2starEchoes,'/SrT2star_uncorr_fin.nii']) 

        fprintf('\n T2 maps already smoothed!..\n')

else

% fprintf('\n Masking ADC with B0 mask...\n')

% before coregistering I smooth the T2map and the T2starmap:
smootfactor=3;
% Smoothing T2map and savings

V = spm_vol([pathT2mapEchoes,'/T2map_fin.nii']);
[T2, ~] = spm_read_vols(V);   

% Normalizzo i valori delle mappe T2 e T2* con massimo 150;
maxT2=max(T2(:));
T2_temp=(T2./maxT2)*150;

% and redefinning threshold
THR=150;

% Salvo
spm_write_vol(V,T2_temp);  

mask1 = (T2>0)>0;
sT2 = smooth3(T2.*mask1,'gaussian',smootfactor);
% mask again after smoothing
sT2 = sT2.*mask1;

    V1.fname = [pathT2mapEchoes,'/ST2map_fin.nii'];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,sT2);  


% Smoothing T2starmap and savings
clear V1

V1 = spm_vol([pathT2starEchoes,'/rT2star_uncorr_fin.nii']);
[T2star, ~] = spm_read_vols(V1);   

% Normalizzo i valori delle mappe T2 e T2* con massimo 150;
maxT2star=max(T2star(:));
T2star_temp=(T2star./maxT2star)*150;

% Salvo
spm_write_vol(V1,T2star_temp);  

mask1 = (T2star>0)>0;
sT2star = smooth3( T2star.*mask1,'gaussian',smootfactor);
% mask again after smoothing
sT2star = sT2star.*mask1;

    V1.fname = [pathT2starEchoes,'/SrT2star_uncorr_fin.nii'];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,sT2star);

end

%% Prima di copiare coregistro le T2map e T2*map alla t1.

if isfile([pathT2mapEchoes,'/T2map_fin_cor_masked.nii.gz']) && isfile([pathT2starEchoes,'/rT2star_uncorr_fin_cor_masked.nii.gz']) 

        fprintf('\n CBV map and CBF map already processed..\n')

else

system(['/usr/local/fsl/bin/flirt -in ',pathT2starEchoes,'/rT2star_uncorr_fin.nii -ref ',path_FLAIR,'/rFLAIR_bet_norm.nii -out ',pathT2starEchoes,'/rT2star_uncorr_fin_cor.nii -omat ',pathT2starEchoes,'/1EchoT2_to_T1.mat -dof 6 -cost mutualinfo -interp trilinear']);
system(['/usr/local/fsl/bin/flirt -in ',pathT2mapEchoes,'/T2map_fin.nii -ref ',path_FLAIR,'/rFLAIR_bet_norm.nii -applyxfm -init ',pathT2starEchoes,'/1EchoT2_to_T1.mat',' -out ',pathT2mapEchoes,'/T2map_fin_cor.nii'])

system(['/usr/local/fsl/bin/fslmaths ' pathT2starEchoes '/rT2star_uncorr_fin_cor.nii -mul ' path_T1w '/T1_Mask.nii ' pathT2starEchoes '/rT2star_uncorr_fin_cor_masked.nii']);
system(['/usr/local/fsl/bin/fslmaths ' pathT2mapEchoes '/T2map_fin_cor.nii -mul ' path_T1w '/T1_Mask.nii ' pathT2mapEchoes '/T2map_fin_cor_masked.nii']);

end

%% Creation of folder with final Coregistered imges

if isfolder([path_general,'OEF/']) 

        fprintf('\n Q-Bold maps already calculated..\n')
        path_OEF=([path_general,'OEF/']);

else

fprintf('\n Creating OEF folder for final calculation...\n')

cd(path_general)
mkdir('OEF')
path_OEF=([path_general,'OEF/']);

copyfile([path_T1w,'/T1_Mask.nii'],[path_general,'OEF/rT1w_bet_mask.nii']);
copyfile([path_FLAIR,'/rFLAIR_bet.nii'],[path_general,'OEF/rrFLAIR_bet.nii']);
copyfile([pathDSCMaps,'/rCBV_masked.nii'],[path_general,'OEF/rrCBV_dsc_new.nii']); % Prendo la CBV corretta
copyfile([pathDSCMaps,'/rCBF_masked.nii'],[path_general,'OEF/rrCBF_dsc_new.nii']); % Prendo la CBV corretta
copyfile([pathT2mapEchoes,'/T2map_fin_cor_masked.nii.gz'],[path_general,'OEF/T2map_fin.nii.gz']);
copyfile([pathT2starEchoes,'/rT2star_uncorr_fin_cor_masked.nii.gz'],[path_general,'OEF/rT2star_uncorr_fin.nii.gz']);
copyfile([pathdiffusion,'/rrADC_masked.nii'],[path_general,'OEF/rrADC.nii']);
system(['gunzip ',path_general,'OEF/*'])
clear list

%% Creation of R2 and R2star from T2 and T2 star and smoothing 3mm R2 e R2 star
% and redefinning threshold
THR=150;
fprintf('\n Calculating R2...\n')

%     mask_V = spm_vol([path_OEF,'/rT1w_bet_mask.nii']);
%     [mask, ~] = spm_read_vols(mask_V); 
%     mask=logical(mask);
    V = spm_vol([path_OEF,'/T2map_fin.nii']);
    [T2, ~] = spm_read_vols(V);    
    T2( isnan(T2)) = 0;
    mask=(T2>0)>0;
    T2 = T2.*mask;
    T2( T2 < 0) = 0;
    T2( T2 > THR) = THR;

    % Save coregistered, cleaned and thresholded maps
    V1.fname=[path_OEF,'/T2map_fin.nii'];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,T2); 
	%----------------------------------------------------------------------
	% Calculate R2
	%----------------------------------------------------------------------
	R2 = 1000 / (T2+eps); %  -> sec^-1
    R2=R2.*mask;
    % Upper boundary to avoid unreasonably bright hotspot
	R2max = 5*median(R2(mask));
    R2(R2>R2max) = R2max;
    R2 = R2.*mask;
	%sR2 = smooth3( R2.*mask1,'gaussian',smootfactor);
	% mask again after smoothing
	%sR2 = R2.*mask1;
    %----------------------------------------------------------------------
    % Save coregistered, cleaned and thresholded maps 
    %----------------------------------------------------------------------
    % R2
    V1.fname = [path_OEF,'R2.nii'];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,R2);  
    %sR2
    %V1.fname = [path_OEF,'sR2.nii'];
%spm_write_vol(V1,sR2); 
    
%% T2*
    %----------------------------------------------------------------------
	% load coregistered T2* map
	%----------------------------------------------------------------------
    V = spm_vol([path_OEF,'/rT2star_uncorr_fin.nii']);
    [T2S, ~] = spm_read_vols(V);
    T2S( isnan( T2S)) = 0;
    mask =(T2S>0)>0;
    T2S = T2S.*mask;
	T2S( T2S < 0) = 0;
	T2S( T2S > THR) = THR;
    V1.fname = [path_OEF,'rT2star_uncorr_fin.nii'];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,T2S); 

    %----------------------------------------------------------------------
	% Calculate R2*
	%----------------------------------------------------------------------
	fprintf('\n Calculating R2S...\n')
   
    R2S = 1000 / (T2S+eps); % -> sec^-1
    % Upper boundary to avoid unreasonably bright hotspot
    R2S = R2S.*mask;
	R2Smax = 5*median(R2S(mask));
	R2S( R2S > R2Smax) = R2Smax;
    R2S = R2S.*mask;
% 	sR2S = smooth3( R2S.*mask,'gaussian',smootfactor);
% 	% mask again after smoothing
% 	sR2S = R2S.*mask;%sR2S = sR2S.*mask;
    
    %---------------------------------------------------------------------
    % Save coregistered, cleaned and thresholded maps 
    %----------------------------------------------------------------------
    % R2S
    V1.fname = [path_OEF,'R2S.nii'];
    spm_write_vol(V1,R2S); 
    %sR2S
%     V1.fname = [path_OEF,'sR2S.nii'];
%     spm_write_vol(V1,sR2S); 
    
	%----------------------------------------------------------------------
	% Calculate smoothed R2'
	%----------------------------------------------------------------------
 	fprintf('\n Calculating R2S-R2 (Rstricht)...\n')
    sR2strich = (R2S - R2);
  	sR2strich(sR2strich < 0) = min(sR2strich(sR2strich>0));
    %sR2strich(isnan( sR2strich)) = 0;
    sR2strich(isnan( sR2strich)) = 0;

    % Upper boundary to avoid unreasonably bright hotspot
	sR2strichmax = 5*median(sR2strich(mask));
%     if sR2strichmax < 20
%     	sR2strichmax = 20; % FIXED threshold
%     end
    sR2strich( sR2strich > sR2strichmax) = sR2strichmax;
	sR2strich = sR2strich.*mask;
    
    %---------------------------------------------------------------------
    % Save coregistered, cleaned and thresholded maps 
    %----------------------------------------------------------------------
    %sR2strich sR2strich.nii
    V1.fname = [path_OEF,'sR2strich.nii'];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,sR2strich); 

%% OEF final calculation:

    mask_V = spm_vol([path_T1w,'/T1_Mask.nii']);
    [mask, ~] = spm_read_vols(mask_V); 
    mask=logical(mask);
    
    %  
fprintf('\n Calculating OEF...\n')

        CBV_hd = spm_vol([path_OEF,'/rrCBV_dsc_new.nii']);
        [CBV_DSC, ~] = spm_read_vols(CBV_hd);    
        C = 4/3*267.61918*pi*0.264*0.42*0.85*3; % 317 Hz
      	CBV_DSC(isnan( CBV_DSC)) = 0;
        CBV_DSC(CBV_DSC<0) = 0;
        CBV_DSC = CBV_DSC.*mask;

        rova_hd.fname=[path_OEF,'rrCBV_dsc_new.nii'];
        rova_hd.dim=CBV_hd.dim;
        rova_hd.mat=CBV_hd.mat;
        rova_hd.n=CBV_hd.n;
        rova_hd.dt=CBV_hd.dt;
        spm_write_vol(rova_hd,CBV_DSC); 
        
        disp(['CBV Values is: ',num2str(mean(CBV_DSC(CBV_DSC>0)))]);

        if mean(CBV_DSC(CBV_DSC>0))>1
            rCBV = CBV_DSC/100.0; %CBV fraction, e.g. CBV_WM = 0.015
            disp(['CBV mean value is (After /100): ',num2str(mean(rCBV(rCBV>0)))]);
        elseif mean(CBV_DSC(CBV_DSC>0))>0.1
            rCBV = CBV_DSC/10.0; %CBV fraction, e.g. CBV_WM = 0.015
            disp(['CBV mean value is (After /10): ',num2str(mean(rCBV(rCBV>0)))]);
        else
            disp(['CBV Values is: ',num2str(mean(rCBV(rCBV>0)))]);
        end

        rCBV( rCBV < 0) = 0;
       	rOEF = (sR2strich ./ (C*rCBV+eps)) .* mask;
        rOEF( rOEF < 0) = 0;
        rOEF(isnan(rOEF))=0;
      	% reasonable boundary to avoid unreasonably high hot spots
        %rOEFmax_DSC = 4*mean(rOEF(rOEF>0 & rOEF<4000));
        rOEF_percentage=rOEF;
        
        Upper_OEF=15;
        rOEFmax_DSC = Upper_OEF;
     	rOEF( rOEF > rOEFmax_DSC) = rOEFmax_DSC;
        rOEF = rOEF .*mask;
    
%         if rOEFmax_DSC < 1.5
%             rOEFmax_DSC = 1.5; % FIXED Threshold
%         end

    rOEFPO2_tmp=rOEF;
    % metto il minimo al posto degli zeri 
    OEF_minimum=min(rOEF(rOEF > 0));
    
    idx_mask=mask==1;
    idx_oef=rOEF==0;

    common_idx1=intersect(idx_mask,idx_oef);
    common_idx=common_idx1==2;

    rOEF(common_idx)=OEF_minimum;
    %rOEF(rOEF==0)=OEF_minimum;
    %rOEF=(rOEF.*mask);%*(10^7);

    %---------------------------------------------------------------------
    % Save coregistered, cleaned and thresholded maps 
    %----------------------------------------------------------------------
    %sR2strich
    CBV_hd.fname =[path_OEF,'OEF.nii'];
    rova_hd.fname=[path_OEF,'OEF.nii'];
    rova_hd.dim=CBV_hd.dim;
    rova_hd.mat=CBV_hd.mat;
    rova_hd.n=CBV_hd.n;
    rova_hd.dt=CBV_hd.dt;
    spm_write_vol(rova_hd,rOEF);   

    % OEF in %
    massimo=max(rOEF(:));
    rOEFperc=(rOEF/massimo)*100;
    rOEFperc(isnan(rOEFperc))=0;
    rova_hd.fname=[path_OEF,'OEF%.nii'];
    rova_hd.dim=CBV_hd.dim;
    rova_hd.mat=CBV_hd.mat;
    rova_hd.n=CBV_hd.n;
    rova_hd.dt=CBV_hd.dt;

    spm_write_vol(rova_hd,rOEFperc); 

%     spm_check_registration([path_OEF,'sR2.nii'],[path_OEF,'sR2S.nii'],[path_OEF,'sR2strich.nii'],[path_OEF,'rrCBV_dsc_new.nii'],[path_OEF,'OEF%.nii']);


%% CMRO2 final calculation:

fprintf('\n Calculating CMRO2...\n')

%OEF*CBF_dsc*C with C=8.68 mmol/ml

rOEF_hd = spm_vol([path_general,'/OEF/OEF.nii']);
[rOEF, ~] = spm_read_vols(rOEF_hd); 

%load CBF
CBF_hd = spm_vol([path_general,'/OEF/rrCBF_dsc_new.nii']);
[CBF, ~] = spm_read_vols(CBF_hd); 
CBF(isnan( CBF)) = 0;
CBF( CBF < 0) = 0;

%Calcolo
CMRO2=(rOEF.*CBF)*8.68;
CMRO2(isnan(CMRO2))=0;
CMRO2(CMRO2 < 0)=0;
% reasonable boundary to avoid unreasonably high hot spots
%      	CMRO2max_DSC = 3*mean(CMRO2(CMRO2>1 & CMRO2<mm));
% %         if rOEFmax_DSC < 1.5
% %             rOEFmax_DSC = 1.5; % FIXED Threshold
% %         end
%      	rOEF( rOEF > rOEFmax_DSC) = rOEFmax_DSC;
%         rOEF = rOEF .*mask;
mm=(max(CMRO2(:)));
CMRO2max_DSC = 5*mean(CMRO2(CMRO2>0));

CMRO2( CMRO2 > CMRO2max_DSC) = CMRO2max_DSC;
CMRO2 = CMRO2 .*mask;

    % metto il minimo al posto degli zeri 
    clear common_idx common_idx1
    CMRO2_minimum=min(CMRO2(CMRO2 > 0));

    idx_CMRO2=CMRO2==0;

    common_idx1=intersect(idx_mask,idx_CMRO2);
    common_idx=common_idx1==2;

    CMRO2(common_idx)=CMRO2_minimum;

%     CMRO2_minimum=min(CMRO2(CMRO2 > 0));
%     CMRO2(CMRO2==0)=CMRO2_minimum;
%     CMRO2=(CMRO2.*mask);%*(10^7);

    rova_hd.fname=[path_OEF,'CMRO2.nii'];
    rova_hd.dim=CBF_hd.dim;
    rova_hd.mat=CBF_hd.mat;
    rova_hd.n=CBF_hd.n;
    rova_hd.dt=CBF_hd.dt;

spm_write_vol(rova_hd,CMRO2);

    massimo=max(CMRO2(:));
    CMRO2perc=(CMRO2/massimo)*100;
    CMRO2perc(isnan(CMRO2perc))=0;
    CMRO2perc(CMRO2perc < 0)=0;
    rova_hd.fname=[path_OEF,'CMRO2%.nii'];
    rova_hd.dim=CBV_hd.dim;
    rova_hd.mat=CBV_hd.mat;
    rova_hd.n=CBV_hd.n;
    rova_hd.dt=CBV_hd.dt;

    spm_write_vol(rova_hd,CMRO2perc); 

%% PO2 final calculation:
% 27(mmHg)*(nthroot(((2./rOEF)-1),2.7))-(CMRO2/4.4)

CC=((2./(rOEFPO2_tmp)));
CC(isnan(CC))=0;
CC(isinf(CC))=0;
CC(CC<0)=0;

CCmax_DSC = 4*mean(CC(CC>0 & CC<300));
CC( CC > CCmax_DSC) = CCmax_DSC;

AA=27*(nthroot(CC,2.7));
AA=(AA/max(AA(:)))*100;
PO2=AA-(CMRO2perc/4.4);
PO2(isnan(PO2))=0;
PO2(PO2<0)=0;
%PO2 = PO2 .*maskrOEF;

% PO2max_DSC = 4*mean(PO2(PO2>0 & PO2<2000));
% PO2( PO2 > PO2max_DSC) = PO2max_DSC;
%PO2 = PO2 .*mask;

    clear common_idx common_idx1
    PO2_minimum=min(PO2(PO2 > 0));

    idx_PO2=PO2==0;

    common_idx1=intersect(idx_mask,idx_PO2);
    common_idx=common_idx1==2;

    PO2(common_idx)=PO2_minimum;

    PO2(isnan(PO2))=0;

CBV_hd.fname = [path_OEF,'PO2.nii'];
    rova_hd.fname=[path_OEF,'PO2.nii'];
    rova_hd.dim=CBV_hd.dim;
    rova_hd.mat=CBV_hd.mat;
    rova_hd.n=CBV_hd.n;
    rova_hd.dt=CBV_hd.dt;
spm_write_vol(rova_hd,PO2);

    massimo=max(PO2(:));
    PO2perc=(PO2/massimo)*100;
    PO2perc(isnan(PO2perc))=0;
    
CBV_hd.fname = [path_OEF,'PO2%.nii'];
    rova_hd.fname=[path_OEF,'PO2%.nii'];
    rova_hd.dim=CBV_hd.dim;
    rova_hd.mat=CBV_hd.mat;
    rova_hd.n=CBV_hd.n;
    rova_hd.dt=CBV_hd.dt;
spm_write_vol(rova_hd,PO2perc);

end

%% Meteonina

if isfolder([path_general,'Clustering/']) 
    path_clustering=[path_general,'Clustering/'];
else
    cd(path_general)
    mkdir('Clustering')
    path_clustering=[path_general,'Clustering/'];
end

if isfolder([path_general,'Metionina/']) 
else
% Coregisternig Metionina
fprintf('\n Coregistering Methionina to standard space!\n')

command_MET_bet=['hd-bet -i ', pathMet,'/Meteonina.nii -device cpu -mode fast -tta 0'];
system(command_MET_bet)

delete([pathMet,'/Meteonina_bet_mask.nii'])
system(['gunzip ',pathMet,'/*'])
movefile([pathMet,'/Meteonina.nii'],[pathMet,'/Not_bet_Meteonina.nii'])
movefile([pathMet,'/Meteonina_bet.nii'],[pathMet,'/Meteonina.nii'])

% Coreg_Met_to_T1(pathMet,path_T1w)

system(['/usr/local/fsl/bin/flirt -in ',pathMet,'/Meteonina.nii -ref ',path_T1w,'/T1w_bet_norm.nii -out ',pathMet,'/rMeteonina.nii']);
system(['gunzip ',pathMet,'/*'])
% Applico maschera T1

command=['/usr/local/fsl/bin/fslmaths ' pathMet '/rMeteonina.nii -mul ' path_T1w '/T1_Mask.nii ' pathMet '/rMeteonina_masked.nii'];
system(command);
system(['gunzip  -f ',pathMet,'/*'])

if isfolder([path_general,'Metionina/']) 
    
    copyfile([pathMet,'/rMeteonina_masked.nii'], [path_general,'/Metionina/rMeteonina_masked.nii'])   
   
else
    mkdir([path_general,'Metionina/'])
    copyfile([pathMet,'/rMeteonina_masked.nii'], [path_general,'/Metionina/rMeteonina_masked.nii'])
end

% Creo metionina masked e la metto nella cartella cluster

%  mask=[patht1segment,'/3DTumor_FLAIR.nii'];
%  im2=[path_general,'/Metionina/rMeteonina_masked.nii'];
% 
%  command=['/usr/local/fsl/bin/fslmaths ' im2 ' -mul ' mask ' ' path_clustering 'Flair_mask_MET.nii'];
%  system(command);
% 
%  system(['gunzip ',path_clustering,'Flair_mask_MET.nii.gz'])

end

%% CLUSTERING 

%Rimpiazzare valori uguali a zero al minimo:

%Moltiplicare per maschera flair:

% Multiply image for FLAIR mask (FSL):

% Adjusting zeros value of Qbold maps:
% zeros values trasformed in minim OEF value and multiplied 10^7;

if isfile([path_clustering,'Flair_mask_ADC.nii']) 

     fprintf('\n Files for clustering already created..\n') 

else

path_OEF=([path_general,'OEF/']);

fprintf('\n Masking ADC, OEF, VP, CMRO2, PO2 and Met with FLAIR mask...\n')

fprintf('\n Placing file within the Clustering folder...\n')
path_clustering=[path_general,'Clustering/'];

%AdjustOEF(path_OEF,patht1segment,path_clustering);

ADC_img=dir([path_OEF,'rrADC.nii']);
OEF_img=dir([path_OEF,'OEF%.nii']);
Vp_img=dir([pathDCE,'/rrVP_masked.nii']);
CMRO_img=dir([path_OEF,'CMRO2%.nii']);
PO2_img=dir([path_OEF,'PO2%.nii']);
Met_img=dir([path_general,'/Metionina/rMeteonina_masked.nii']);

tomul={[ADC_img.folder,'/',ADC_img.name], [OEF_img.folder,'/',OEF_img.name],[Vp_img.folder,'/',Vp_img.name],[CMRO_img.folder,'/',CMRO_img.name],[PO2_img.folder,'/',PO2_img.name],[Met_img.folder,'/',Met_img.name]};
outname={'Flair_mask_ADC.nii','Flair_mask_OEF.nii','Flair_mask_Vp.nii','Flair_mask_CMRO.nii','Flair_mask_PO2.nii','Flair_mask_MET.nii'};

%% MOLTIPLICA PER LA MAPPA rOEF La maschera in modo che esca solo la regione compresa

V = spm_vol([path_OEF,'/rT2star_uncorr_fin.nii']);
[T2S, ~] = spm_read_vols(V);
ss=size(T2S);
T2S(isnan(T2S))=0;
WW_tmp=(T2S>0)>0;

disp('Here you need the 3DTumor_FLAIR_reor (whole tumor segmentation)!')


if isfile([patht1segment,'3DTumor_FLAIR_reor.nii']) 
    
    %system(['gunzip ',patht1segment,'/3DTumor_FLAIR_reor.nii.'])
    V2= spm_vol([patht1segment,'/3DTumor_FLAIR_reor.nii']);
    [tumormask, ~] = spm_read_vols(V2);
    %system(['gzip ',patht1segment,'/3DTumor_FLAIR_reor.nii.gz'])

else
    
    V2= spm_vol([patht1segment,'/3DTumor_FLAIR.nii']);
    [tumormask, ~] = spm_read_vols(V2);

end

mask_intersection=WW_tmp.*tumormask;

for i=1:length(tomul)

    V = spm_vol(tomul{i});
    [Current, ~] = spm_read_vols(V);
   
    Current=Current.*mask_intersection;
    
    Current(isnan(Current))=0;

    V1.fname = [path_clustering,'/',outname{i}];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,Current);

    clear current_minimum Current V V2 Current_WT

end
 
    fprintf('\n FLAIR Segmentation is absent!\n')

end


end
fprintf('\n *****************************\n')
fprintf('\n *********** DONE ************\n')
fprintf('\n *****************************\n')
fprintf('\n Now run the Clustering Analysis! \n')
fprintf('\n Cheers\n')
fprintf('\n *****************************\n')
fprintf('\n *****************************\n')
fprintf('\n *****************************\n')