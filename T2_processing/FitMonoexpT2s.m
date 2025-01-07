function [outputArg1,outputArg2] = FitMonoexpT2s(T2map_TR,T2star_TR,TET2,TET2s,pathT2mapEchoes,pathT2starEchoes,THR)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
%% Fit T2s:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%% Fit T2 parameter MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g2_val=4000;
thr=THR;
fprintf('\n Starting T2s Processing:\n');
fprintf('\n T2 map Processing:\n');

    cd(pathT2mapEchoes);
    list=dir('T*nii');
    nechoes = length(list); % assume that number of files in list equals number of echoes
    %Read first data set from list and determine parameters
    headerinfo_T2 = ...
    spm_vol([pathT2mapEchoes,'/',list(1).name]); 
    [vol, ~] = spm_read_vols( headerinfo_T2);
    [nx, ny, nslices] = size(vol);
    clear vol;

    % read T2 echo datat and parameters
    vol(1:nx,1:ny,1:nslices,1:nechoes) = 0;
    TE(1:nechoes) = 0;
    TR(1:nechoes) = 0;
    % read json file for TE & TR info
    for i = 1:nechoes
        headerinfo_T2 = ...
        spm_vol([pathT2mapEchoes,'/',list(i).name]); 
        newStr_tmp = extractAfter(headerinfo_T2.descrip,"=");
        %TE(i) = str2double(extractBefore(newStr_tmp,";")); %[in ms]
        TR(i) = T2map_TR;  % [in ms]
        headerinfo_T2 = ...
            spm_vol([pathT2mapEchoes,'/',list(i).name]); 
        [vol(:,:,:,i), ~] = spm_read_vols( headerinfo_T2);
        %vol(:,:,:,i) = smooth3(vol(:,:,:,i),'gaussian',3);
    end

    % ---------------------------------------------------------------------
    % Save nifti Header as template for saving of Maps (TE1, T2, errors):
    % ---------------------------------------------------------------------
    V_save       = headerinfo_T2(1); % save image header info as structure
    V_save.dt    = [64 0];           % 64 stands for float data type 
    
    %----------------------------------------------------------------------
    % Calculate T2 maps
    %----------------------------------------------------------------------
    %maske = squeeze(vol(:,:,:,1)) > mean(mean(mean(vol(:,:,:,1))));
        fprintf('\n3D T2 data: Fitting all echoes!\n');
        zeit = TET2;
        data = vol;
        clear TE vol;
    
    fprintf('\n Starting T2 fit:\n slice# ');

    % slice wise fit of T2 and amp values
    for i = 1:nslices
        fprintf('%d ', i);
        serie = squeeze(data(:,:,i,:));
        % function 
        %[amap,cmap,devmap] = ...
        %       imagefit2param(zeit,serie,x,y,g1,g2,p,noisefakt,maxdev)
        % 2 parameter fit:  wert[i] = a*{x+y*exp(-zeit[i]/c)}
        % zeit: vector containing time points
        % serie: 3D set containing weighted images
        x = 0;          % x according to fit formula
        y = 1;          % y according to fit formula
        g1 = 1;         % lower limit for c
        g2 = g2_val;      % upper limit for c
        p = 0.01;       % accuracy in relative values
        noisefakt = 1; % noise factor
        maxdev = 20;    % max. deviation of data from fitted values [in %]
        [~, cmap, devmap] = ...
            imagefit2param( zeit, serie,x,y,g1,g2,p,noisefakt,maxdev);
        t2map(:,:,i) = cmap;
        error(:,:,i) = devmap;
    end

    fprintf('\n');

    % mask T2 map and restrict maximum value to 150 ms to improve results
    % of coregistration / spatial normalisation
    %t2map = t2map.*maske;
    t2map(t2map>thr)=thr;
    %error = error.*maske;
    
    % save T2 parameter map
    fprintf('\n Saving T2 maps\n');

    V_save.fname = [pathT2mapEchoes,'/T2map_fin.nii'];
    spm_write_vol(V_save, t2map);

    % save T2 fit error map
     V_save.fname = [pathT2mapEchoes,'/T2map_error_fin.nii'];
    spm_write_vol(V_save, error);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%% Fit T2* parameter MAP %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear list nx ny nslices pathT2mapEchoes data zeit
    nechoes = 10;
    cd(pathT2starEchoes);
    list=dir('T*nii');
    
    %Read first data set from list and determine parameters
    headerinfo_T2 = ...
    spm_vol([pathT2starEchoes,'/',list(1).name]); 
    [vol, ~] = spm_read_vols( headerinfo_T2);
    [nx, ny, nslices, ~]  = size(vol);
    clear vol;
    
    % read json file for TE & TR info
    for i = 1:nechoes
        headerinfo_T2 = ...
        spm_vol([pathT2starEchoes,'/',list(i).name]); 
        newStr_tmp = extractAfter(headerinfo_T2.descrip,"=");
        %TE(i) = str2double(extractBefore(newStr_tmp,";")); %[in ms]
        TR(i) = T2star_TR;  % [in ms]
        headerinfo_T2 = ...
            spm_vol([pathT2starEchoes,'/',list(i).name]); 
        [vol(:,:,:,i), ~] = spm_read_vols( headerinfo_T2);
        %vol(:,:,:,i) = smooth3(vol(:,:,:,i),'gaussian',3);
    end
    
    V_save       = headerinfo_T2(1); % save image header info as structure
    V_save.dt    = [64 0];           % 64 stands for float data type 
    
    data=vol;
    zeit = TET2s;
	T2S_uncorr(1:nx,1:ny,1:nslices) = 0;
	T2S_error(1:nx,1:ny,1:nslices) = 0;
    fprintf('\n3D T2 star data: Fitting all echoes!\n');
    
    fprintf('\n Starting T2 fit:\n slice# ');
    % slice wise fit of T2 and amp values
    for i = 1:nslices
        fprintf('%d ', i);
        serie = squeeze(data(:,:,i,:));
        % function 
        %[amap,cmap,devmap] = ...
        %       imagefit2param(zeit,serie,x,y,g1,g2,p,noisefakt,maxdev)
        % 2 parameter fit:  wert[i] = a*{x+y*exp(-zeit[i]/c)}
        % zeit: vector containing time points
        % serie: 3D set containing weighted images
        x = 0;          % x according to fit formula
        y = 1;          % y according to fit formula
        g1 = 1;         % lower limit for c
        g2 = g2_val;      % upper limit for c
        p = 0.01;       % accuracy in relative values
        noisefakt = 10; % noise factor
        maxdev = 20;    % max. deviation of data from fitted values [in %]
        [~, cmap, devmap] = ...
            imagefit2param( zeit, serie,x,y,g1,g2,p,noisefakt,maxdev);
        T2S_uncorr(:,:,i) = cmap;
        T2S_error(:,:,i) = devmap;
    end

    fprintf('\n');

    % mask T2 map and restrict maximum value to 150 ms to improve results
    % of coregistration / spatial normalisation
     %maske = squeeze(data(:,:,:,1)) > mean(mean(mean(data(:,:,:,1))));
     %T2S_uncorr = T2S_uncorr.*maske;
     T2S_uncorr(T2S_uncorr>thr)=thr;
     %T2S_error = T2S_error.*maske;
    
    % save uncorrected T2* parameter map
    fprintf('\n Saving T2 star maps\n');

        V_save.fname = [pathT2starEchoes,'/T2star_uncorr_fin.nii'];
        spm_write_vol(V_save, T2S_uncorr);

        V_save.fname = [pathT2starEchoes,'/T2star_uncorr_err_fin.nii'];
        spm_write_vol(V_save, T2S_uncorr);
end

