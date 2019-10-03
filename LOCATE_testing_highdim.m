function LOCATE_testing_highdim(test_image_directory_name, varargin)
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

addpath('MATLAB');

training_image_directory_name = pwd;

if nargin > 0
    training_image_directory_name = varargin{1};
end

% Assigning the Root directories
root_data_directory = test_image_directory_name;
results_directory = sprintf('%s/LOCATE_results_directory',test_image_directory_name);

% Creating the results directory
if ~exist(results_directory)
    mkdir(results_directory)
end

xdir = dir(sprintf('%s/*_BIANCA_LPM.nii.gz',train_image_directory_name));
if numel(xdir) == 0
    error('Cannot find any input image. Please check your training_image_directoy_name ');
else
    xdir = dir(sprintf('%s/*_BIANCA_LPM.nii.gz',root_data_directory));
    xsplit = regexp(xdir(1).name,'_BIANCA_LPM','split');
    xfeats = dir(sprintf('%s/%s_feature_*',root_data_directory,xsplit{1}));
    if numel(xfeats) == 0
        fprintf('Warning: No feature modality found. Can extract only volume and ventricle distance \n');
    end
end

numfeats = numel(xfeats);

feature_selection_cols = ones(numfeats+2,1);
if nargin > 2
    if numel(varargin{2}) == 1
        error('Third input (feature_select) must be a vector with number of elements equal to the number of features specified.');
    elseif numel(varargin{2}) < numfeats + 2
        error('Number of columns in feature_select does not match the number of features specified.');
    else
        feature_selection_cols = varargin{2};
    end
end

verbose = 0;
if nargin > 3
    verbose = varargin{3};
end

if verbose
    feature_selection_cols
    train_image_directory_name
end

% Loading the Trained Regression model
try
    load(sprintf('%s/LOCATE_training_files/RF_regression_model_LOCATE.mat',train_image_directory_name));
catch
    error('Training model not found. Run LOCATE_training first or check the training images directory');
end

if sum(feature_selection_cols) ~= (size(RFmodel.X,2)/20)
    error('The number of features selected does not match the features used to train the model. Kindly check your feature_select input or use a different model.');
end

xdir = dir(sprintf('%s/*_BIANCA_LPM.nii.gz',root_data_directory));
imgfeatmats = cell(numel(xdir),1);
ventdistfeatmats = cell(numel(xdir),1);
lesvolfeatmats = cell(numel(xdir),1);
index_indices_list = cell(numel(xdir),1);
index_maps = cell(numel(xdir),1);

for subj = 1:numel(xdir)    
    
    if verbose
        subj
    end
    
    xsplit = regexp(xdir(subj).name,'_BIANCA_LPM','split');
    xfeats = dir(sprintf('%s/%s_feature_*',root_data_directory, xsplit{1}));
    
    does_file_exist = exist(sprintf('%s/LOCATE_features_%s.mat',results_directory,xsplit{1}),'file');
    
    if does_file_exist == 2
        load(sprintf('%s/LOCATE_features_%s.mat',results_directory,xsplit{1}))
    else
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
        [flairintfeats, ventdistfeats, lesvolfeats, index_numbers, index_mask] ...
            = LOCATE_test_feature_extraction(lesionmask, ventdistmap, flairimage, index_mask, index_numbers);
        if verbose
            fprintf('LOCATE features extracted! \n')
        end
    end
    
    flairintfeats = single(flairintfeats);
    ventdistfeats = single(ventdistfeats);
    lesvolfeats = single(lesvolfeats);
    index_numbers = single(index_numbers);
    index_mask = single(index_mask);
    
    % Storing the features in a cell array
    imgfeatmats{subj} = flairintfeats;
    ventdistfeatmats{subj} = ventdistfeats;
    lesvolfeatmats{subj} = lesvolfeats;
    
    index_indices_list{subj} = index_numbers;
	index_maps{subj} = index_mask;
    
    % Saving the features as an intermediate results .mat file
    save(sprintf('%s/LOCATE_features_%s.mat',...
        results_directory,xsplit{1}),'flairintfeats',...
    'ventdistfeats','lesvolfeats','minbestthr_values','index_mask');
    
end
% % Saving the array as an intermediate results .mat file
save(sprintf('%s/LOCATE_features.mat',results_directory),'imgfeatmats',...
    'ventdistfeatmats','lesvolfeatmats','index_indices_list','index_maps');

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

% Initializing cell array to store images
% lesion_masks = cell(1,numel(xdir));
% thresholds = cell(1,numel(xdir));
% thresholdmap = cell(1,numel(xdir));
% indexmaps = cell(1,numel(xdir));

for testsubj = 1:numel(xdir)
    % Load the test image
    xsplit = regexp(xdir(testsubj).name,'_BIANCA_LPM','split');
    lesionmaskfile = sprintf('%s/%s_BIANCA_LPM.nii.gz',root_data_directory,xsplit{1}); 
    lesionmask = read_avw(lesionmaskfile);
    
    % Extract the corresponding test features and index masks from cell
    % arrays
    testflairfeatmat = imgfeatmats{testsubj};
    testventdistfeatmat = ventdistfeatmats{testsubj};
    testlesvolfeatmat = lesvolfeatmats{testsubj};
    index_mask = index_maps{testsubj};
    index_numbers = index_indices_list{testsubj};
    voronoi_test_features_all = [testventdistfeatmat,testlesvolfeatmat,testflairfeatmat];
    if verbose
        fprintf('Voronoi test features has been collected! \n');
    end
    threshold_array = 0:0.05:0.95;
    feature_selection_cols_exp = repmat(feature_selection_cols, [numel(threshold_array),1]);
    feature_selection_cols_exp = feature_selection_cols_exp(:)';
    voronoi_test_features = voronoi_test_features_all(:,feature_selection_cols_exp>0);
    testmeanbestthrs = predict(RFmodel,voronoi_test_features);
    
    %Assigning the values to the final image
    final_binary_lesionmask = zeros(size(lesionmask));
    threshold_mask = zeros(size(lesionmask));
    
    for ind = 1:numel(index_numbers)  
        voronoi_region_mask = index_mask == index_numbers(ind);
        lesionmask_voronoi = lesionmask.*double(voronoi_region_mask);
        binary_voronoi_lesionmask = lesionmask_voronoi > testmeanbestthrs(ind);
        thresh = single(voronoi_region_mask).*testmeanbestthrs(ind);
        threshold_mask = threshold_mask + thresh;
        final_binary_lesionmask = final_binary_lesionmask | binary_voronoi_lesionmask;
    end
    lesion_masks{testsubj} = final_binary_lesionmask;
    thresholds{testsubj} = testmeanbestthrs;
    thresholdmap{testsubj} = threshold_mask;
    indexmaps{testsubj} = index_mask;
    
    % Saving the images
    save_avw(index_mask,sprintf('%s/%s_indexmap.nii.gz',results_directory,xsplit{1}),'f',[1 1 1]);
    save_avw(threshold_mask,sprintf('%s/%s_thresholdsmap.nii.gz',results_directory,xsplit{1}),'f',[1 1 1]);
    save_avw(final_binary_lesionmask,sprintf('%s/%s_BIANCA_LOCATE_binarylesionmap.nii.gz',results_directory,xsplit{1}),'f',[1 1 1]);
    save(sprintf('%s/%s_LOCATE_thresholds.mat',results_directory,xsplit{1}),'testmeanbestthrs');
    
    copying_image_geometry(lesionmaskfile, sprintf('%s/%s_indexmap.nii.gz',results_directory,xsplit{1}), verbose);
    copying_image_geometry(lesionmaskfile, sprintf('%s/%s_thresholdsmap.nii.gz',results_directory,xsplit{1}), verbose);
    copying_image_geometry(lesionmaskfile, sprintf('%s/%s_BIANCA_LOCATE_binarylesionmap.nii.gz',results_directory,xsplit{1}), verbose);
    
end
% Save all the results as a consolidated single folder
save(sprintf('%s/Consolidated_LOCATE_output.mat',results_directory),'lesion_masks','thresholds','indexmaps','thresholdmap');

