function LOCATE_testing_per_subject(LPM_filename, training_model_filename, features_filename, feature_select, num_of_modalities, results_directory)
% Function for testing LOCATE
%   Copyright - FMRIB, WIN, University of Oxford
%   Vaanathi Sundaresan - 25/05/2018
% 
%   Mandatory input:
%   test_image_directory_name - directory consisting of test features
%
%   Example funtional calls:
%   1. LOCATE_testing(test_image_directory_name, train_image_directory_name);
%    - Name of the directory where the test and training images are located
%   2. LOCATE_testing(test_image_directory_name, train_image_directory_name, feature_select);
%    - If you want to select specific features for training and testing
%   3. LOCATE_testing(test_image_directory_name, train_image_directory_name, feature_select, verbose);


addpath('MATLAB');

if numel(feature_select) ~= num_of_modalities+2
    error('Please check the number of elements in feature_select');
end
try
    load(training_model_filename);
catch
    error('Training model not found. Run LOCATE_training first or check the training images directory');
end

if sum(feature_select) ~= (size(RFmodel.X,2)/20)
    error('The number of features selected does not match the features used to train the model. Kindly check your feature_select input or use a different model.');
end

% Creating the results directory
if ~exist(results_directory,'dir')
    mkdir(results_directory)
end


%if numel(xdir) == 0
%    error('Cannot find any feature file. Please check your test_image_directory_name');
%end

imgfeatmats = cell(1,1);
ventdistfeatmats = cell(1,1);
lesvolfeatmats = cell(1,1);

xsplit = regexp(features_filename,'LOCATE_features_','split');
does_file_exist = exist(sprintf('%s/%s_indexmap.nii.gz',results_directory,xsplit{2}(1:end-4)),'file');

if does_file_exist ~= 2
    % Load the test image
    lesionmask = read_avw(LPM_filename);

    % Load the test feature
    matfile = load(features_filename);
    flairintfeats = matfile.flairintfeats;
    ventdistfeats = matfile.ventdistfeats;
    lesvolfeats = matfile.lesvolfeats;
    index_numbers = matfile.index_numbers;
    index_mask = matfile.index_mask;

    imgfeatmats{1} = flairintfeats;
    ventdistfeatmats{1} = ventdistfeats;
    lesvolfeatmats{1} = lesvolfeats;


    [flairfeatmats_modified,ventdistfeatmats_modified,lesvolfeatmats_modified] = LOCATE_modify_voronoi_features(imgfeatmats,ventdistfeatmats,lesvolfeatmats);
    imgfeatmats = flairfeatmats_modified;
    ventdistfeatmats = ventdistfeatmats_modified;
    lesvolfeatmats = lesvolfeatmats_modified;

    testflairfeatmat = imgfeatmats{1};
    testventdistfeatmat = ventdistfeatmats{1};
    testlesvolfeatmat = lesvolfeatmats{1};

    voronoi_test_features_all = [testventdistfeatmat,testlesvolfeatmat,testflairfeatmat];
    fprintf('Voronoi test features has been collected! \n');

    threshold_array = 0:0.05:0.95;
    feature_selection_cols_exp = repmat(feature_select, [numel(threshold_array),1]);
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

    % Saving the images
    save_avw(index_mask,sprintf('%s/%s_indexmap.nii.gz',results_directory,xsplit{2}(1:end-4)),'f',[1 1 1]);
    save_avw(threshold_mask,sprintf('%s/%s_thresholdsmap.nii.gz',results_directory,xsplit{2}(1:end-4)),'f',[1 1 1]);
    save_avw(final_binary_lesionmask,sprintf('%s/%s_BIANCA_LOCATE_binarylesionmap.nii.gz',results_directory,xsplit{2}(1:end-4)),'f',[1 1 1]);
    save(sprintf('%s/%s_LOCATE_thresholds.mat',results_directory,xsplit{2}(1:end-4)),'testmeanbestthrs');

    verbose = 1;
    copying_image_geometry(LPM_filename, sprintf('%s/%s_indexmap.nii.gz',results_directory,xsplit{2}(1:end-4)), verbose);
    copying_image_geometry(LPM_filename, sprintf('%s/%s_thresholdsmap.nii.gz',results_directory,xsplit{2}(1:end-4)), verbose);
    copying_image_geometry(LPM_filename, sprintf('%s/%s_BIANCA_LOCATE_binarylesionmap.nii.gz',results_directory,xsplit{2}(1:end-4)), verbose);
end      

end
