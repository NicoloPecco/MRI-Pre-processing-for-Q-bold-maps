function  WhiteStrip_Norm(path_T1w,path_FLAIR,patht1segment)


%% T1 NORM

% Loading tumor region:

V = spm_vol([patht1segment,'/3DTumor_FLAIR.nii']);
[Tumor, ~] = spm_read_vols(V);
Tumor(isnan(Tumor))=0;

% Loading anat mask: 
V = spm_vol([path_T1w,'/T1_Mask.nii']);
[Mask, ~] = spm_read_vols(V);
Mask(isnan(Mask))=0;

% Loading T1 image:

VT = spm_vol([path_T1w,'/T1w_bet.nii']);
[T1w, ~] = spm_read_vols(VT);
T1w(isnan(T1w))=0;

% Loading T1 White matter segmentation:
V = spm_vol([path_T1w,'/c2T1w_orig.nii']);
%V = spm_vol([path_T1w,'/c2T1w.nii']);
[WM, ~] = spm_read_vols(V);
WM(isnan(WM))=0;

% Calculating Mean and std:
WM_Bin=(WM>0.9)>0;
WM_S=WM_Bin-Tumor;

t1_mul=(T1w.*WM_S);
M_t1=mean(t1_mul(t1_mul>0));
S_t1=std(t1_mul(t1_mul>0));

T1w_new = (T1w - M_t1)/S_t1;
T1w_new=(T1w_new.*Mask);
T1w_new=T1w_new+abs(min(T1w_new(:)));
T1w_new=(T1w_new.*Mask);

%% FLAIR NORM

% Loading FLAIR image:

VF = spm_vol([path_FLAIR,'/rFLAIR_bet.nii']);
[Flair, ~] = spm_read_vols(VF);
Flair(isnan(Flair))=0;

Flair_mul=(Flair.*WM_S);
M_t1=mean(Flair_mul(Flair_mul>0));
S_t1=std(Flair_mul(Flair_mul>0));

Flair_new = (Flair - M_t1)/S_t1;
Flair_new=(Flair_new.*Mask);
Flair_new=Flair_new+abs(min(Flair_new(:)));
Flair_new=(Flair_new.*Mask);

% Salvo anatomiche normalizzate 

VT.fname=[path_T1w,'/T1w_bet_norm.nii'];
VT.dt(1)=16;
spm_write_vol(VT,T1w_new);

VF.fname=[path_FLAIR,'/rFLAIR_bet_norm.nii'];
VF.dt(1)=16;
spm_write_vol(VF,Flair_new);

return
end