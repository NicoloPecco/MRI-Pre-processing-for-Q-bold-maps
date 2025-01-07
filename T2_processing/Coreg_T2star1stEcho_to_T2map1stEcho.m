function [outputArg1,outputArg2] = Coreg_T2star1stEcho_to_T2map1stEcho(pathT2mapEchoes,pathT2starEchoes)
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
    
    % Taking T2 map ( GLOBAL REFERENCE)
	PG = spm_select('FPList', [pathT2mapEchoes], '^T2map_fin.nii$');
    if isempty(PG)
        error('No Target (T2, TE1) selected. Check files or file names!');
    end      
% 	PG = PG(1,:);
	VG = spm_vol(PG);

    % SOURCE: T2 star map
	PF = spm_select( 'FPList', [pathT2starEchoes], '^T2star_uncorr_f.*.nii$');  
	if isempty(PF)
     	error('No Source (coregistered T2star) selected. Check files or file names!');
	end
    PF = PF(1,:);
 	VF = spm_vol(PF);

    % select "other images" for coregistration (FLAIR segments)
	PO1 = spm_select( 'FPList', [pathT2starEchoes], '^T2star_uncorr_e.*.nii$');
    if isempty(PO1)
    	error('No tissue segments selected. Check files or file names!');
    end
    
	if isempty(PO1) && isempty(PO2) && isempty(PO3)
      	PO = PF;
	else
    	PO = char(PF,PO1);
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

    V = spm_vol([pathT2starEchoes,'/rT2star_uncorr_fin.nii']);
    [T2S, ~] = spm_read_vols(V);
    T2S( isnan( T2S)) = 0;
    V1.fname = [pathT2starEchoes,'/rT2star_uncorr_fin.nii'];
    V1.dim=V.dim;
    V1.mat=V.mat;
    V1.n=V.n;
    V1.dt=V.dt;
    spm_write_vol(V1,T2S); 

end