function [outputArg1,outputArg2] = Coreg_Rstrich_to_T1(path_FLAIR,pathT2mapEchoes,pathT2starEchoes,path_OEF)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here
	global defaults;
	spm_get_defaults;
	flags = defaults.coreg;
    	% Flags to pass to routine (spm_reslice) to create resliced images
    resFlags = struct(...
    	'interp', 1,...%flags.write.interp,... % 1: trilinear interpolation
    	'wrap', flags.write.wrap,...           % wrapping info (ignore...)
      	'mask', flags.write.mask,...           % masking (see spm_reslice)
      	'which',1,...                          % write reslice time series for 
      	'mean',0);                             % later use do not write 
    
    % Taking T1
%  	PG = spm_select('FPList',[path_FLAIR] , '^rFLAIR.*bet.nii$'); 
    PG = spm_select('FPList',[path_FLAIR] , '^rFLAIR.*.nii$');
    if isempty(PG)
        error('No Target T1w selected. Check files or file names!');
    end      
	PG = PG(1,:);
	VG = spm_vol(PG);

    % SOURCE: First echo t2map
	PF = spm_select('FPList', [pathT2mapEchoes], '^T2map.*.nii$');
	if isempty(PF)
     	error('No Source T2map selected. Check files or file names!');
	end
    PF = PF(10,:);
 	VF = spm_vol(PF);

    %select "other images" for coregistration 
	PO1 = spm_select( 'FPList', [path_OEF], '^sR.*.nii$');
    if isempty(PO1)
    	error('No sR2strich selected. Check files or file names!');
    end
    
  	PO2 = spm_select( 'FPList', [path_OEF], '^R.*.nii$');
    if isempty(PO2)
    	error('No sR2strich selected. Check files or file names!');
    end

    PO3 = spm_select( 'FPList', [path_OEF], '^rT2.*.nii$');
    if isempty(PO3)
    	error('No sR2strich selected. Check files or file names!');
    end

    PO4 = spm_select( 'FPList', [path_OEF], '^T2m.*.nii$');
    if isempty(PO4)
    	error('No sR2strich selected. Check files or file names!');
    end

	if isempty(PO1) && isempty(PO2) 
      	PO = PF;
    else
        PO = char(PF,PO1,PO2,PO3,PO4);
	end

	% do coregistration, i.e. get coregistration parameters
	x  = spm_coreg( VG, VF, flags.estimate);

	% get the transformation to be applied with parameters 'x'
	M  = inv( spm_matrix(x));

	% in MM we put the transformations for the 'other' images
	MM = zeros( 4, 4, size( PO,1));

	for j=1:size( PO,1)
        % get the transformation matrix for every image
      	MM(:,:,j) = spm_get_space( deblank(PO(j,:)));
    end

	for j=1:size( PO,1)
        % write the transformation applied to every image
      	% filename: deblank(PO(j,:))
      	spm_get_space(deblank(PO(j,:)), M*MM(:,:,j));
	end

 	P = char( PG, PO);
 	spm_reslice( P,resFlags);
end

