function [flag] = process_anatomy(path_T1w,path_FLAIR,path_general)
%UNTITLED3 Summary of this function goes here
%   Processing of anatomical data (T1w and FLAIR)
listT1=dir(path_T1w);
listT1(ismember({listT1.name},{'.','..','.DS_Store'}))=[];

if length(listT1)>=11
fprintf('\n T1-w Processing already done..\n');
else
fprintf('\n Starting T1-w Processing \n');
PF = spm_select('ExtFPList',[path_T1w],'.*g.nii$');
PF = PF(1,:);
segment_SPM12(PF, 6);
fprintf('\n T1-w Processing done!\n');
%fprintf('\n Starting FLAIR Processing');
end

% listFL=dir(path_FLAIR);
% listFL(ismember({listFL.name},{'.','..','.DS_Store'}))=[];
% 
% if length(listFL)>=11
% fprintf('\n FLAIR Processing already done..\n');
% else
% PF = spm_select('ExtFPList',[path_FLAIR],'.*FLAIR.nii$');
% PF = PF(1,:);
% segment_SPM12(PF, 6);
% fprintf('\n FLAIR Processing done..\n');
% flag=1;
% end

% if isfile([path_T1w,'/T1w_bet.nii']);
%     fprintf('\n T1w already betted!\n');
% else
% fprintf('\n Creating T1 brain mask...\n')
% 
% command=['bet2 ',path_T1w,'/T1w.nii ',path_T1w, '/T1w_bet.nii -m -f 0.5'];
% system(command);
% system(['gunzip ',path_T1w,'/*'])
% end
% 
% if isfile([path_FLAIR,'/FLAIR_bet.nii'])
%         fprintf('\n FLAIR already betted!\n');
% else
% fprintf('\n Creating FLAIR brain mask...\n')
% 
% command=['bet2 ',path_FLAIR,'/FLAIR.nii ',path_FLAIR, '/FLAIR_bet.nii -m'];
% system(command);
% system(['gunzip ',path_FLAIR,'/*'])
% end


%TumorLesionSegmentation(path_FLAIR,path_general); % Ottengo TumorLesion.nii

% Coregistering TumorLesion.nii to T1w.nii for visualization porpuse:

% global defaults;
% spm_get_defaults;
% flags = defaults.coreg;
% flags.write.interp=0;
% % Flags to pass to routine (spm_reslice) to create resliced images
% resFlags = struct(...
% 'interp', 0,...%flags.write.interp,... % 1: trilinear interpolation
% 'wrap', flags.write.wrap,...           % wrapping info (ignore...)
% 'mask', flags.write.mask,...           % masking (see spm_reslice)
% 'which',1,...                          % write reslice time series for
% 'mean',0);                             % later use do not write
% % Taking T1
% PG = spm_select('FPList', [path_T1w], '^T1w.nii$');
% if isempty(PG)
% error('No Target (T2, TE1) selected. Check files or file names!');
% end
% % 	PG = PG(1,:);
% VG = spm_vol(PG);
% % SOURCE: T2 star map
% PF = spm_select( 'FPList', [path_FLAIR], '^FLAIR.*.nii$');
% if isempty(PF)
% error('No Source (coregistered T2star) selected. Check files or file names!');
% end
% PF = PF(1,:);
% VF = spm_vol(PF);
% % select "other images" for coregistration (FLAIR segments)
% copyfile([path_general,'/TumorLesion/TumorLesion.nii'],[path_general,'/TumorLesion/T1_TumorLesion.nii']);
% PO1 = spm_select( 'FPList', [path_general,'/TumorLesion'], '^T1_TumorLe.*.nii$');
% if isempty(PO1)
% error('No tissue segments selected. Check files or file names!');
% end
% 
% PO = char(PO1);
% 
% x  = spm_coreg( VG, VF, flags.estimate);
% % get the transformation to be applied with parameters 'x'
% M  = inv( spm_matrix(x));
% % in MM we put the transformations for the 'other' images
% MM = zeros( 4, 4, size( PO,1));
% for j=1:size( PO,1)
% % get the transformation matrix for every image
% MM(:,:,j) = spm_get_space( deblank(PO(j,:)));
% end
% 
% for j=1:size( PO,1)
% % write the transformation applied to every image
% % filename: deblank(PO(j,:))
% spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
% end
% 
% P = char( PG, PO);
% spm_reslice( P,resFlags);
% 
% % Creating unique segmentation image for visualization:
% 
%         V = spm_vol([path_T1w,'/c1T1w.nii']);
%         [GM, ~] = spm_read_vols(V);
%     
%         % FLAIR WM    
%         V = spm_vol([path_T1w,'/c2T1w.nii']);
%         [WM, ~] = spm_read_vols(V);
%         
%         % FLAIR CSF     
%         V = spm_vol([path_T1w,'/c3T1w.nii']);
%         [CSF, ~] = spm_read_vols(V);
%         
%         % FLAIR Tumor  
%         V = spm_vol([path_general,'/TumorLesion/rT1_Tumorlesion.nii']);
%         [TR, ~] = spm_read_vols(V);
%         
%         totseg=CSF+GM*2+WM*3+TR*4;
%         V.fname =[path_general,'/TumorLesion/All_seg_T1.nii'];
%         spm_write_vol(V,totseg); 
        
end

