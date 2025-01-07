function [TET2,TET2s,flagT2,flagT2s] = Motion_Coregistr_T2s(pathT2mapEchoes,pathT2starEchoes,motT2)

% Motion correction of echoes in T2 and T2 star indipendentely

% Coregistration of echoes and risclice of T2star on T2  (reduce B0 inhom. effects).

global defaults;
	spm_get_defaults;
	flags = defaults.coreg;
    	% Flags to pass to routine (spm_reslice) to create resliced images
    resFlags = struct(...
    	'interp', 0,...%flags.write.interp,... % 0: NN
    	'wrap', flags.write.wrap,...           % wrapping info (ignore...)
      	'mask', flags.write.mask,...           % masking (see spm_reslice)
      	'which',1,...                          % write reslice time series for 
      	'mean',0);          
%% T2
prefix='T2map';
T2_filepath = spm_select( 'FPList', pathT2mapEchoes, ...
['^' prefix '.*\.nii$']);

if motT2(1)==1
mkdir([pathT2mapEchoes,'/Orig_vol']); 
for i=1:size(T2_filepath,1)
    headerinfo_T2=spm_vol([T2_filepath(i,:)]); 
    newStr_tmp = extractAfter(headerinfo_T2.descrip,"=");
    TET2(i) = str2double(extractBefore(newStr_tmp,";")); 
    movefile([T2_filepath(i,:)],[pathT2mapEchoes,'/Orig_vol/'])
end
% Run Realignment: choose 'estwrite', i.e. estimate and write
VG = spm_select( 'FPList', [pathT2mapEchoes,'/Orig_vol'], ...
'^T2map_e1.*.nii$');

for jj=1:size(T2_filepath,1)-1

VF= spm_select( 'FPList', [pathT2mapEchoes,'/Orig_vol'], ...
['^T2map_e',num2str(jj+1),'.*.nii$']);

% do coregistration, i.e. get coregistration parameters
	x  = spm_coreg( VG, VF, flags.estimate);

	% get the transformation to be applied with parameters 'x'
	M  = inv( spm_matrix(x));

	% in MM we put the transformations for the 'other' images
	MM = zeros( 4, 4);

    % get the transformation matrix for every image
    MM= spm_get_space( deblank(VF(1,:)));

    % write the transformation applied to every image
    % filename: deblank(PO(j,:))
    spm_get_space(deblank(VF), M*MM);
 	P = char( VG, VF);
    
 	spm_reslice(P,resFlags);
end

    for i=1:size(T2_filepath,1)
    if i==1
    copyfile([pathT2mapEchoes,'/Orig_vol/T2map_e',num2str(i),'.nii'],[pathT2mapEchoes,'/T2map_e',num2str(i),'.nii'])
    else
    movefile([pathT2mapEchoes,'/Orig_vol/rT2map_e',num2str(i),'.nii'],[pathT2mapEchoes,'/T2map_e',num2str(i),'.nii'])
    end
    end
    
    fprintf('\nMotion correction of T2 echoes performed\n')
    flagT2=1;
    
else
    for i=1:size(T2_filepath,1)
    headerinfo_T2=spm_vol([T2_filepath(i,:)]); 
    newStr_tmp = extractAfter(headerinfo_T2.descrip,"=");
    TET2(i) = str2double(extractBefore(newStr_tmp,";")); 
    end
    fprintf('\nMotion correction of T2 echoes not performed\n')
    flagT2=0;
end

%% T2s
clear T2_filepath headerinfo_T2 VG VF M MM P
prefix='T2star';
T2_filepath = spm_select( 'FPList', pathT2starEchoes, ...
['^' prefix '.*\.nii$']);
if motT2(2)==1
mkdir([pathT2starEchoes,'/Orig_vol']);
for i=1:size(T2_filepath,1)
    headerinfo_T2=spm_vol([T2_filepath(i,:)]); 
    newStr_tmp = extractAfter(headerinfo_T2.descrip,"=");
    TET2s(i) = str2double(extractBefore(newStr_tmp,";")); 
    movefile([T2_filepath(i,:)],[pathT2starEchoes,'/Orig_vol/'])
end
% Run Realignment: choose 'estwrite', i.e. estimate and write
T2_filepath = spm_select( 'FPList', [pathT2starEchoes,'/Orig_vol'], ...
['^' prefix '.*\.nii$']);
% Run Realignment: choose 'estwrite', i.e. estimate and write
VG = spm_select( 'FPList', [pathT2starEchoes,'/Orig_vol'], ...
'^T2star_e01.*.nii$');

for jj=1:size(T2_filepath,1)-1

VF= spm_select( 'FPList', [pathT2starEchoes,'/Orig_vol'], ...
['^T2star_e',num2str(jj+1,'%02d'),'.*.nii$']);

% do coregistration, i.e. get coregistration parameters
	x  = spm_coreg( VG, VF, flags.estimate);

	% get the transformation to be applied with parameters 'x'
	M  = inv( spm_matrix(x));

	% in MM we put the transformations for the 'other' images
	MM = zeros( 4, 4);

    % get the transformation matrix for every image
    MM= spm_get_space( deblank(VF(1,:)));

    % write the transformation applied to every image
    % filename: deblank(PO(j,:))
    spm_get_space(deblank(VF), M*MM);
 	P = char( VG, VF);
    
 	spm_reslice( P,resFlags);
end 

% sposto i file motion corrected nella cartella principale senza prefisso

    for i=1:size(T2_filepath,1)
    if i==1
    copyfile([pathT2starEchoes,'/Orig_vol/T2star_e',num2str(i,'%02d'),'.nii'],[pathT2starEchoes,'/T2star_e',num2str(i,'%02d'),'.nii'])
    else
    movefile([pathT2starEchoes,'/Orig_vol/rT2star_e',num2str(i,'%02d'),'.nii'],[pathT2starEchoes,'/T2star_e',num2str(i,'%02d'),'.nii'])
    end
    end
    fprintf('\nMotion correction of T2 echoes performed\n')
    flagT2s=1;
%      resFlags = struct(...
%     	'interp', 1,...%flags.write.interp,... % 1: Trilinear
%     	'wrap', flags.write.wrap,...           % wrapping info (ignore...)
%       	'mask', flags.write.mask,...           % masking (see spm_reslice)
%       	'which',1,...                          % write reslice time series for 
%       	'mean',0);     
%         % Coregistration: T2star echoes to T2map first echo in order to:
% 
%                     % Work in the same space:
%                     % Reduce the effect of B0 inhogeneities (Trilinear interp)
% 
%         VG = spm_select( 'FPList', [pathT2mapEchoes],'^T2map_e1.*.nii$');
% 
%         VF= spm_select( 'FPList', [pathT2starEchoes],['^T2star_e','.*.nii$']);
%         VF=VF(1,:);
% 
%         % select "other images" for coregistration GM, WM, CSF segments & masks
%         PO = spm_select( 'FPList', [pathT2starEchoes], ['^T2star_e','.*.nii$']);
% 
%         % do coregistration
% 
%         % this method from spm_coreg_ui.m
%         % get coregistration parameters
%         x  = spm_coreg( VG, VF, flags.estimate);
% 
%         % get the transformation to be applied with parameters 'x'
%         M  = inv( spm_matrix(x));
% 
%         % in MM we put the transformations for the 'other' images
%         MM = zeros( 4, 4, size( PO,1));
% 
%         for j=1:size( PO,1)
%             % get the transformation matrix for every image
%             MM(:,:,j) = spm_get_space( deblank( PO(j,:)));
%         end
% 
%         for j=1:size( PO,1)
%             % write the transformation applied to every image
%             % filename: deblank(PO(j,:))
%             spm_get_space( deblank( PO(j,:)), M*MM(:,:,j));
%         end
% 
%         P = char(VG, PO);
%         spm_reslice( P,resFlags);
% 
% T2_filepath = spm_select( 'FPList', [pathT2starEchoes,'/Orig_vol'],['^' prefix '.*\.nii$']);        
%         
%     for i=1:size(T2_filepath,1)
% 
%         movefile([pathT2starEchoes,'/T2star_e',num2str(i,'%02d'),'.nii'],[pathT2starEchoes,'/Orig_vol/mot_corr_T2star_e',num2str(i,'%02d'),'.nii'])
%         movefile([pathT2starEchoes,'/rT2star_e',num2str(i,'%02d'),'.nii'],[pathT2starEchoes,'/T2star_e',num2str(i,'%02d'),'.nii'])
% 
%     end

else
    for i=1:size(T2_filepath,1)
    headerinfo_T2=spm_vol([T2_filepath(i,:)]); 
    newStr_tmp = extractAfter(headerinfo_T2.descrip,"=");
    TET2s(i) = str2double(extractBefore(newStr_tmp,";")); 
    end
    fprintf('\nMotion correction of T2 star echoes not performed\n')
    flagT2s=0;
end


end
