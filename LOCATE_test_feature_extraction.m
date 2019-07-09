function [flairintfeats, ventdistfeats, lesvolfeats, index_numbers, index_mask] ...
        = LOCATE_test_feature_extraction(lesionmask, ventdistmap, flairimage, index_mask, index_numbers)
% Function for extracting test features for LOCATE
%   Copyright - FMRIB, WIN, University of Oxford
%   Vaanathi Sundaresan - 25/05/2018 
%
%   lesionmask - BIANCA lesion probability maps
%   ventdistmap - distance from the ventricles
%   flairimage - FLAIR image
%   index_mask - Voronoi tessellation index mask
%   index_numbers - Voronoi tessellation region indices

    thresholds = 0:0.05:0.95;    
    fprintf('Extra details');
    size(lesionmask)
    numel(flairimage)
    size(index_mask)
    numel(thresholds)
    flairintfeats = zeros(numel(index_numbers),numel(thresholds)*numel(flairimage));
    ventdistfeats = zeros(numel(index_numbers),numel(thresholds));
    lesvolfeats = zeros(numel(index_numbers),numel(thresholds));
    % Iterating over Voronoi regions
    for ind = 1:numel(index_numbers)  
        voronoi_region_mask = index_mask == index_numbers(ind);
        lesionmask_voronoi = lesionmask.*double(voronoi_region_mask);
        flairint_values = zeros(numel(flairimage),numel(thresholds));
        ventdist_values = zeros(1,numel(thresholds));
        lesvol_values = zeros(1,numel(thresholds));
        % Applying thresholds and calculating features for each thresholds
        for thr = 1:numel(thresholds)
            binary_lesion_mask = lesionmask_voronoi > thresholds(thr);
            for img_no = 1:numel(flairimage)
                flair_image = flairimage{img_no};
                size(flair_image)
                size(flairint_values)
		size(binary_lesion_mask)
                flairint_values(img_no,thr) = mean(flair_image(binary_lesion_mask));
            end
            ventdist_values(thr) = mean(ventdistmap(binary_lesion_mask));
            lesvol_values(thr) = sum(binary_lesion_mask(:));            
        end
        flairint_values = flairint_values';
        flairint_values = flairint_values(:)';
        % Replacing NANs with -1
        flairint_values(isnan(flairint_values)) = -1;
        ventdist_values(isnan(ventdist_values)) = -1;
        lesvol_values(isnan(lesvol_values)) = -1;
      
        % Aggregating features
        flairintfeats(ind,:) = flairint_values;
        ventdistfeats(ind,:) = ventdist_values;
        lesvolfeats(ind,:) = lesvol_values;
    end
