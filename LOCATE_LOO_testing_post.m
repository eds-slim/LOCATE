function LOCATE_LOO_testing(varargin)
% Function for extracting  features from the images in a directory
% and performing Leave-one-subject-out testing for LOCATE
%   Copyright - FMRIB, WIN, University of Oxford
%   Vaanathi Sundaresan - 25/05/2018
%
%   Example funtional calls:
%   1. LOCATE_LOO_testing(train_image_directory_name);
%    - Name of the directory where images for LOO evaluation are located
%   2. LOCATE_LOO_testing(train_image_directory_name, feature_select);
%    - If you want to select specific features for training and testing
%   3. LOCATE_LOO_testing(train_image_directory_name, feature_select, verbose);
%
%   Optional inputs (in the order):
%    - feature_select - vector with elements indicating if the feature has to be included or not. Current order is distance from ventricles, lesion volume and other modalities in alphabetical naming order
%                  (e.g. If FLAIR is the only modality provided and distance from ventricles is not needed then feature_select = [0, 1, 1])
%    - verbose (0 or 1)


if nargin > 0
    training_image_directory_name = varargin{1};
end

results_directory = sprintf('%s/LOCATE_LOO_results_directory',training_image_directory_name);

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

if verbose
    feature_selection_cols
    training_image_directory_name
end



% Creating the results directory
if ~exist(results_directory, 'dir')
    mkdir(results_directory)
end

imgfeatmats = cell(numel(xdir),1);
ventdistfeatmats = cell(numel(xdir),1);
lesvolfeatmats = cell(numel(xdir),1);
minbestthrs = cell(numel(xdir),1);
maxbestthrs = cell(numel(xdir),1);
meanbestthrs = cell(numel(xdir),1);
index_indices_list = cell(numel(xdir),1);
index_maps = cell(numel(xdir),1);



features_all_fn = sprintf('%s/LOCATE_features.mat',results_directory);

if exist(features_all_fn,'file')
    load(features_all_fn)
else
    
    for subj = 1:numel(xdir)
        
        if verbose
            subj
        end
        
        xsplit = regexp(xdir(subj).name,'_BIANCA_LPM','split');
        features_fn = sprintf('%s/LOCATE_features_%s.mat',results_directory,xsplit{1});
        
        if exist(features_fn,'file')
            load(features_fn)
        else
            error('Extracted features not available. Run *_sub.m first!')
        end
        % Storing the features in a cell array
        imgfeatmats{subj} = flairintfeats;
        ventdistfeatmats{subj} = ventdistfeats;
        lesvolfeatmats{subj} = lesvolfeats;
        minbestthrs{subj} = minbestthr_values;
        maxbestthrs{subj} = maxbestthr_values;
        meanbestthrs{subj} = meanbestthr_values;
        index_indices_list{subj} = index_numbers;
        index_maps{subj} = index_mask;
        
    end
    % Saving the array as an intermediate results .mat file
    save(features_all_fn, 'imgfeatmats',...
        'ventdistfeatmats','lesvolfeatmats','minbestthrs','maxbestthrs','meanbestthrs','index_indices_list','index_maps', '-v7.3');
    
end

if verbose
    fprintf('LOCATE features saved! \n')
end

