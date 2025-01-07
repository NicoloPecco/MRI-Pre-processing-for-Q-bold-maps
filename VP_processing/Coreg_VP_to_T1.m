function [outputArg1,outputArg2] = Coreg_VP_to_T1(pathDCE,path_T1w)
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
	PG = spm_select('FPList', [path_T1w], '^T.*bet.nii$');
    if isempty(PG)
        error('No Target (T2, TE1) selected. Check files or file names!');
    end      
	PG = PG(1,:);
	VG = spm_vol(PG);

    % SOURCE: FLAIR 
	PF = spm_select('FPList',[pathDCE,'/'] , '^B.*bet.nii$');   
	if isempty(PF)
     	error('No Source B0 bet of DCE!');
	end
    PF = PF(1,:);
 	VF = spm_vol(PF);


    PO1 = spm_select( 'FPList', [pathDCE,'/'], '^rVP.*d.nii$');
    if isempty(PO1)
    	disp('No VP selected.');
    end
    
    if isempty(PO1)
      	PO = PF;
    else
        PO = char(PF,PO1);%PO1
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