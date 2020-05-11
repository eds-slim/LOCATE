function LOCATE_training_sub(varargin)
% Function for extracting  features from the images in a directory and training for LOCATE
%   Copyright - FMRIB, WIN, University of Oxford
%   Vaanathi Sundaresan - 25/05/2018
%
%   Example funtional calls:
%   1. LOCATE_training(train_image_directory_name);
%    - Name of the directory where the training images for feature extraction are located
%   2. LOCATE_training(train_image_directory_name, feature_select);
%    - If you want to select specific features for training and testing
%   3. LOCATE_training(train_image_directory_name, feature_select, verbose);
%
%   Optional inputs (in the order):
%    - feature_select - vector with elements indicating if the feature has to be included or not. Current order is distance from ventricles, lesion volume and other modalities in alphabetical naming order
%                 (e.g. If FLAIR is the only modality provided and distance from ventricles is not needed then feature_select = [0, 1, 1])
%    - verbose (0 or 1)

if nargin > 0
    training_image_directory_name = varargin{1};
end

% Assigning the Root directories
results_directory = sprintf('%s/LOCATE_training_files',training_image_directory_name);

% Creating the results directory
if ~exist(results_directory, 'dir')
    mkdir(results_directory)
end

xdir = dir(sprintf('%s/*_BIANCA_LPM.nii.gz',training_image_directory_name));
if numel(xdir) == 0
    error('Cannot find any input image. Please check your training_image_directoy_name');
else
    xsplit = regexp(xdir(1).name,'_BIANCA_LPM','split');
    xfeats = dir(sprintf('%s/%s_feature_*',training_image_directory_name,xsplit{1}));
    if numel(xfeats) == 0
        fprintf('Warning: No feature modality found. Can extract only volume and ventricle distance \n');
    end
end

numfeats = numel(xfeats);

feature_selection_cols = ones(numfeats+2,1);
if nargin > 1
    if numel(varargin{2}) == 1
        error('Second input (feature_select) must be a vector with number of elements equal to the number of features specified.');
    elseif numel(varargin{2}) < numfeats + 2
        error('Number of columns in feature_select does not match the number of features specified.');
    else
        feature_selection_cols = varargin{2};
    end
end

verbose = 0;
if nargin > 2
    verbose = varargin{3};
end

subj = nan;
if nargin > 3
    subj = varargin{4};
end
if isnan(subj); error('subj not set'); end

if verbose
    feature_selection_cols
    training_image_directory_name
end

imgfeatmats = cell(numel(xdir),1);
ventdistfeatmats = cell(numel(xdir),1);
lesvolfeatmats = cell(numel(xdir),1);
minbestthrs = cell(numel(xdir),1);
maxbestthrs = cell(numel(xdir),1);
meanbestthrs = cell(numel(xdir),1);
index_indices_list = cell(numel(xdir),1);
index_maps = cell(numel(xdir),1);

if verbose
    subj
end

xsplit = regexp(xdir(subj).name,'_BIANCA_LPM','split');

features_fn = sprintf('%s/LOCATE_features_%s.mat',results_directory,xsplit{1});
if exist(features_fn,'file')
    load(features_fn)
else
    xfeats = dir(sprintf('%s/%s_feature_*',training_image_directory_name,xsplit{1}));
    flairimage = cell(numel(xfeats),1);
    % Loading the image files
    lesionmaskfile = sprintf('%s/%s_BIANCA_LPM.nii.gz',training_image_directory_name,xsplit{1});
    manualmaskfile = sprintf('%s/%s_manualmask.nii.gz',training_image_directory_name,xsplit{1});
    biancamaskfile = sprintf('%s/%s_biancamask.nii.gz',training_image_directory_name,xsplit{1});
    brainmaskfile = sprintf('%s/%s_brainmask.nii.gz',training_image_directory_name,xsplit{1});
    lesionmask = read_avw(lesionmaskfile);
    manualmask = read_avw(manualmaskfile);
    biancamask = read_avw(biancamaskfile);
    brainmask = read_avw(brainmaskfile);
    
    if feature_selection_cols(1) == 0
        try
            ventdistmapfile = sprintf('%s/%s_ventdistmap.nii.gz',training_image_directory_name,xsplit{1});
            ventdistmap = read_avw(ventdistmapfile);
        catch
            ventdistmap = zeros(size(lesionmask));
        end
    else
        ventdistmapfile = sprintf('%s/%s_ventdistmap.nii.gz',training_image_directory_name,xsplit{1});
        ventdistmap = read_avw(ventdistmapfile);
    end
    
    for subj_feat_no = 1:numel(xfeats)
        flairimagefile = sprintf('%s/%s',training_image_directory_name,xfeats(subj_feat_no).name);
        flairimage{subj_feat_no} = read_avw(flairimagefile);
    end
    
    if verbose
        fprintf('All specified feature image modalities loaded \n');
    end
    
    % Getting image dimensions and determining up/downsampling factor
    dim = size(lesionmask);
    factor = 150./dim; %floor(max(dim)./dim); %
    inv_factor = 1./factor;
    
    % Up/downsampling the images
    lesionmask_resized = imresizen(lesionmask,factor);
    biancamask = imresizen(single(biancamask),factor);
    brainmask = imresizen(single(brainmask),factor);
    biancamask = (biancamask>0) & (brainmask>0);
    
    
    % Performing Voronoi tessellation on resampled images
    [~, index_mask, index_numbers] = LOCATE_Voronoi_tessellation_resized(lesionmask_resized, biancamask, inv_factor);
    if verbose
        fprintf('Voronoi Tessellation done! \n')
    end
    numel(index_numbers)
    index_mask = imresizen(index_mask,inv_factor,'nearest');
    
    if verbose
        fprintf('Resizing of Index mask done! \n')
    end
    index_numbers = setdiff(union(index_mask(:),[]),0);
    numel(index_numbers)
    
    % Extractng features from Voronoi regions individually
    [flairintfeats, ventdistfeats, lesvolfeats, minbestthr_values, maxbestthr_values, meanbestthr_values, index_numbers, index_mask] ...
        = LOCATE_feature_extraction(lesionmask, ventdistmap, flairimage, manualmask, index_mask, index_numbers);
    if verbose
        fprintf('LOCATE features extracted! \n')
    end
    
    
    flairintfeats = single(flairintfeats);
    ventdistfeats = single(ventdistfeats);
    lesvolfeats = single(lesvolfeats);
    minbestthr_values = single(minbestthr_values);
    maxbestthr_values = single(maxbestthr_values);
    meanbestthr_values = single(meanbestthr_values);
    index_numbers = single(index_numbers);
    index_mask = single(index_mask);
    
    % Saving the features as an intermediate results .mat file
    save(features_fn,...
        'flairintfeats', 'ventdistfeats', 'lesvolfeats',...
        'minbestthr_values', 'maxbestthr_values', 'meanbestthr_values',...
        'index_numbers','index_mask');
end


