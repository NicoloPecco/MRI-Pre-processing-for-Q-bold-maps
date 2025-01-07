function Risizeimage(path_T1w,isotropicVoxel)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

if length(dir([path_T1w,'/*']))<7

V = spm_vol([path_T1w,'/T1w_orig.nii']);
[T1, ~] = spm_read_vols(V); 
vx_vol = sqrt(sum(V.mat(1:3,1:3).^2));

        if round(sum(vx_vol),4)==3 
        
        disp('T1w_orig image is alredy isotropic 1x1x1.')
        %movefile([path_T1w,'/T1w_orig.nii'],[path_T1w,'/T1w.nii'])
        
        else
        
        command=['flirt -in ' path_T1w '/T1w_orig.nii  -ref ' path_T1w, '/T1w_orig.nii -out ',path_T1w,'/T1w -applyisoxfm ' num2str(isotropicVoxel)];
        system(command)
        command2=['gunzip ' path_T1w '/*'];
        system(command2)
        delete([path_T1w,'/T1w_orig.nii'])
        movefile([path_T1w,'/T1w.nii'],[path_T1w,'/T1w_orig.nii'])
        end

else

    fprintf('\n T1 Processing already done! \n')

end
