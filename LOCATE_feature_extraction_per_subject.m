function LOCATE_feature_extraction_per_subject(feature_select, num_of_modalities, features_path, subject_name, varargin)

% Function for extracting features for LOCATE
%   Copyright - FMRIB, WIN, University of Oxford
%   Vaanathi Sundaresan - 25/05/2018
% 
%   Mandatory inputs:
%   feature_select (binary integer array) - vector with elements indicating if the feature has to be included or not. Current order is distance from ventricles, lesion volume and other modalities in alphabetical naming order
%                 (e.g. If FLAIR is the only modality provided and distance from ventricles is not needed then feature_select = [0, 1, 1])
%   num_of_modalities (integer) - number of input modalities used (e.g.
%   num_of_modalities = 1, if only FLAIR is used, num_of_modalities = 2 if
%   FLAIR and T1 are used)
%   features_path (string) - path to the directory where features need to be stored
%   subject_name (string) - To save the LOCATE features file
%
%   Example funtional calls:
%   1. LOCATE_testing(feature_select, num_of_modalities, path_to_file1, ..., path_to_fileN );
%   path_to_files - starting from the paths to modalities, LPM from BIANCA,
%   biancamask, brainmask, ventricle distancemap (if the first value of
%   feature_select is 1)
  
addpath('MATLAB');

if ~exist(features_path)
    mkdir(features_path)
end
nargin
if feature_select(1) == 0
    if nargin < (2+num_of_modalities+3)
        error('Number of input paths provided is lesser than those required. Please check the input filepaths provided');
    end
else
    if nargin < (2+num_of_modalities+4)
        error('Number of input paths provided is lesser than those required. Please check the input filepaths provided');
    end
end

flairimage = cell(num_of_modalities,1);
c = 1;
for subj_feat_no = 1:num_of_modalities
    flairimage{subj_feat_no} = read_avw(varargin{c});
    c = c+1;
end

lesionmask = read_avw(varargin{num_of_modalities+1});   
biancamask = read_avw(varargin{num_of_modalities+2});
brainmask = read_avw(varargin{num_of_modalities+3});

if feature_select(1) == 0
    try 
        ventdistmap = read_avw(varargin{num_of_modalities+4});
    catch
        ventdistmap = zeros(size(lesionmask));
    end
else
    ventdistmap = read_avw(varargin{num_of_modalities+4});
end

fprintf('All specified feature image modalities loaded \n');
        
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
fprintf('Voronoi Tessellation done! \n')
index_mask = imresizen(index_mask,inv_factor,'nearest');
index_numbers = setdiff(union(index_mask(:),[]),0);
numel(index_numbers)
% Extractng features from Voronoi regions individually 
[flairintfeats, ventdistfeats, lesvolfeats, index_numbers, index_mask] ...
    = LOCATE_test_feature_extraction(lesionmask, ventdistmap, flairimage, index_mask, index_numbers);
        
fprintf('LOCATE features extracted! \n')        
    
flairintfeats = single(flairintfeats);
ventdistfeats = single(ventdistfeats);
lesvolfeats = single(lesvolfeats);
%minbestthr_values = single(minbestthr_values);
%maxbestthr_values = single(maxbestthr_values);
%meanbestthr_values = single(meanbestthr_values);
index_numbers = single(index_numbers);
index_mask = single(index_mask);

% Saving the features as an intermediate results .mat file
save(sprintf('%s/LOCATE_features_%s.mat',features_path,subject_name),'flairintfeats',...
'ventdistfeats','lesvolfeats','index_numbers','index_mask');
    
end
