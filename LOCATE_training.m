function LOCATE_training(varargin)
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

addpath('MATLAB');

training_image_directory_name = pwd;

if nargin > 0
    training_image_directory_name = varargin{1};
end
% Assigning the Root directories
root_data_directory = training_image_directory_name;
results_directory = sprintf('%s/LOCATE_training_files',training_image_directory_name);

% Creating the results directory
if ~exist(results_directory)
    mkdir(results_directory)
end

xdir = dir(sprintf('%s/*_BIANCA_LPM.nii.gz',root_data_directory));
if numel(xdir) == 0
    error('Cannot find any input image. Please check your training_image_directoy_name');
else
    xsplit = regexp(xdir(1).name,'_BIANCA_LPM','split');
    xfeats = dir(sprintf('%s/%s_feature_*',root_data_directory,xsplit{1}));
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
    
imgfeatmats = cell(numel(xdir),1);
ventdistfeatmats = cell(numel(xdir),1);
lesvolfeatmats = cell(numel(xdir),1);
minbestthrs = cell(numel(xdir),1);
maxbestthrs = cell(numel(xdir),1);
meanbestthrs = cell(numel(xdir),1);
index_indices_list = cell(numel(xdir),1);
index_maps = cell(numel(xdir),1);

for subj = 1:numel(xdir)
    if verbose
        subj
    end
    
    xsplit = regexp(xdir(subj).name,'_BIANCA_LPM','split');
    
    does_file_exist = exist(sprintf('%s/LOCATE_features_%s.mat',results_directory,xsplit{1}),'file');
    
    if does_file_exist == 2
        load(sprintf('%s/LOCATE_features_%s.mat',results_directory,xsplit{1}))
    else
        xfeats = dir(sprintf('%s/%s_feature_*',root_data_directory,xsplit{1}));
        flairimage = cell(numel(xfeats),1);
        % Loading the image files
        lesionmaskfile = sprintf('%s/%s_BIANCA_LPM.nii.gz',root_data_directory,xsplit{1});
        manualmaskfile = sprintf('%s/%s_manualmask.nii.gz',root_data_directory,xsplit{1});    
        biancamaskfile = sprintf('%s/%s_biancamask.nii.gz',root_data_directory,xsplit{1});
        brainmaskfile = sprintf('%s/%s_brainmask.nii.gz',root_data_directory,xsplit{1});
        lesionmask = read_avw(lesionmaskfile);
        manualmask = read_avw(manualmaskfile);    
        biancamask = read_avw(biancamaskfile);
        brainmask = read_avw(brainmaskfile);
        
        if feature_selection_cols(1) == 0   
            try 
                ventdistmapfile = sprintf('%s/%s_ventdistmap.nii.gz',root_data_directory,xsplit{1});
                ventdistmap = read_avw(ventdistmapfile);
            catch
                ventdistmap = zeros(size(lesionmask));
            end
        else
            ventdistmapfile = sprintf('%s/%s_ventdistmap.nii.gz',root_data_directory,xsplit{1});
            ventdistmap = read_avw(ventdistmapfile);
        end
        
        for subj_feat_no = 1:numel(xfeats)
            flairimagefile = sprintf('%s/%s',root_data_directory,xfeats(subj_feat_no).name);
            flairimage{subj_feat_no} = read_avw(flairimagefile);
        end 
        
        if verbose
            fprintf('All specified feature image modalities loaded \n');
        end
        
        % Getting image dimensions and determining up/downsampling factor
        dim = size(lesionmask);
        factor = 150./dim;%floor(max(dim)./dim);
        inv_factor = 1./factor;
        
        % Up/downsampling the images
        lesionmask = imresizen(lesionmask,factor);
        biancamask = imresizen(single(biancamask),factor);
        brainmask = imresizen(single(brainmask),factor);
        biancamask = (biancamask>0) & (brainmask>0);  
    
        % Performing Voronoi tessellation on resampled images
        [lesionmask, index_mask, index_numbers] = LOCATE_Voronoi_tessellation_resized(lesionmask, biancamask, inv_factor);
        if verbose
            fprintf('Voronoi Tessellation done! \n')
        end
        numel(index_numbers)
        % Extractng features from Voronoi regions individually 
        [flairintfeats, ventdistfeats, lesvolfeats, minbestthr_values, maxbestthr_values, meanbestthr_values, index_numbers, index_mask] ...
            = LOCATE_feature_extraction(lesionmask, ventdistmap, flairimage, manualmask, index_mask, index_numbers);
        if verbose
            fprintf('LOCATE features extracted! \n')
        end
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
    
    % Saving the features as an intermediate results .mat file
    save(sprintf('%s/LOCATE_features_%s.mat',...
        results_directory,xsplit{1}),'flairintfeats',...
    'ventdistfeats','lesvolfeats','minbestthr_values','maxbestthr_values','meanbestthr_values','index_numbers','index_mask');
    
end
% Saving the array as an intermediate results .mat file
save(sprintf('%s/LOCATE_features.mat',results_directory),'imgfeatmats',...
    'ventdistfeatmats','lesvolfeatmats','minbestthrs','maxbestthrs','meanbestthrs','index_indices_list','index_maps');

if verbose
    fprintf('LOCATE features saved! \n')
end
% Loading the features (not needed) if run on a local system (can aslo
% resumed from this step onwards if you already have features saved in above format
% as 'LOCATE_features.mat'
load(sprintf('%s/LOCATE_features.mat',results_directory))

% Modifying features to correct it to remove -1's (internal processing)
[flairfeatmats_modified,ventdistfeatmats_modified,lesvolfeatmats_modified] = LOCATE_modify_voronoi_features(imgfeatmats,ventdistfeatmats,lesvolfeatmats);
imgfeatmats = flairfeatmats_modified;
ventdistfeatmats = ventdistfeatmats_modified;
lesvolfeatmats = lesvolfeatmats_modified;
if verbose
    fprintf('Voronoi features has been corrected! \n');
end
% Concatenating the training features
trainflairfeatmat = [];
trainventdistfeatmat = [];
trainlesvolfeatmat = [];
trainminbestthrs = [];
trainmaxbestthrs = [];
trainmeanbestthrs = [];
for trainsubj = 1:numel(xdir)
    trainflairfeatmat = [trainflairfeatmat;imgfeatmats{trainsubj}];
    trainventdistfeatmat = [trainventdistfeatmat;ventdistfeatmats{trainsubj}];
    trainlesvolfeatmat = [trainlesvolfeatmat;lesvolfeatmats{trainsubj}];
    trainminbestthrs = [trainminbestthrs;minbestthrs{trainsubj}];
    trainmaxbestthrs = [trainmaxbestthrs;maxbestthrs{trainsubj}];
    trainmeanbestthrs = [trainmeanbestthrs;meanbestthrs{trainsubj}]; 
end
voronoi_train_features_all = [trainventdistfeatmat,trainlesvolfeatmat,trainflairfeatmat];
if verbose
    fprintf('Voronoi training features has been collected! \n');
end
threshold_array = 0:0.05:0.95;
feature_selection_cols_exp = repmat(feature_selection_cols, [numel(threshold_array),1]);
feature_selection_cols_exp = feature_selection_cols_exp(:)';
voronoi_train_features = voronoi_train_features_all(:,feature_selection_cols_exp>0);

% Training RF regression model
RFmodel = TreeBagger(1000,voronoi_train_features,trainmaxbestthrs,'Method','Regression',...
    'numPredictorsToSample','all');
if verbose
    fprintf('Training done! \n');
end
% Saving the model
save(sprintf('%s/RF_regression_model_LOCATE.mat',results_directory),'RFmodel');
if verbose
    fprintf('LOCATE model saved! \n');
end