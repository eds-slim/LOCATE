function [flairintfeats, ventdistfeats, lesvolfeats, minbestthr_values, maxbestthr_values, meanbestthr_values, index_numbers, index_mask] ...
        = LOCATE_feature_extraction(lesionmask, ventdistmap, flairimage, manualmask, index_mask, index_numbers)
% Function for extracting test features for LOCATE
%   Copyright - FMRIB, WIN, University of Oxford
%   Vaanathi Sundaresan - 25/05/2018
%
%   lesionmask - BIANCA lesion probability maps
%   ventdistmap - distance from the ventricles
%   flairimage - FLAIR image
%   index_mask - Voronoi tessellation index mask
%   manualmask - sumal segmentation
%   index_numbers - Voronoi tessellation region indices        
    thresholds = 0:0.05:0.95;
    
    flairintfeats = zeros(numel(index_numbers),numel(thresholds)*numel(flairimage));
    ventdistfeats = zeros(numel(index_numbers),numel(thresholds));
    lesvolfeats = zeros(numel(index_numbers),numel(thresholds));
    dice_values = zeros(numel(index_numbers),1);
    maxbestthr_values = zeros(numel(index_numbers),1);
    minbestthr_values = zeros(numel(index_numbers),1);
    meanbestthr_values = zeros(numel(index_numbers),1);
    fprintf('Initiating LOCATE feature extraction...')
    % Iterating over Voronoi regions
    for ind = 1:numel(index_numbers)
        voronoi_region_mask = index_mask == index_numbers(ind);
        lesionmask_voronoi = lesionmask.*double(voronoi_region_mask);
        manualmask_voronoi = single(manualmask).*double(voronoi_region_mask);
        flairint_values = zeros(numel(flairimage),numel(thresholds));
        ventdist_values = zeros(1,numel(thresholds));
        lesvol_values = zeros(1,numel(thresholds));
        size(lesvol_values)
        diceindex_values = zeros(1,numel(thresholds));
        size(flairint_values)
        % Applying thresholds and calculating features and dice indices for each thresholds
        for thr = 1:numel(thresholds)
            binary_lesion_mask = lesionmask_voronoi > thresholds(thr);
            size(binary_lesion_mask)
            size(flairimage)
            for img_no = 1:numel(flairimage)
                flair_image = flairimage{img_no};
                flairint_values(img_no,thr) = mean(flair_image(binary_lesion_mask));
            end
            ventdist_values(thr) = mean(ventdistmap(binary_lesion_mask));
            lesvol_values(thr) = sum(binary_lesion_mask(:));
            intersect_img = binary_lesion_mask & manualmask_voronoi;
            no_intersect = sum(intersect_img(:));
            diceindex_values(thr) = 2*(no_intersect/(sum(manualmask_voronoi(:))+sum(binary_lesion_mask(:))));
            if sum(manualmask_voronoi(:)) == 0 && sum(binary_lesion_mask(:))== 0
                diceindex_values(thr) = 1;
            end
        end
        flairint_values = flairint_values';
        flairint_values = flairint_values(:)';
        % Replacing NANs with -1
        diceindex_values(isnan(diceindex_values)) = 0;
        flairint_values(isnan(flairint_values)) = -1;
        ventdist_values(isnan(ventdist_values)) = -1;
        lesvol_values(isnan(lesvol_values)) = -1;
        [maxdice] = max(diceindex_values);
        maxdicepos = diceindex_values == max(diceindex_values);
        dice_values(ind) = maxdice;
        % Aggregating features and (mean, min and max) thresholds that give maximum dices
        maxbestthr_values(ind) = max(thresholds(maxdicepos));
        minbestthr_values(ind) = min(thresholds(maxdicepos));
        meanbestthr_values(ind) = mean(thresholds(maxdicepos));
        flairintfeats(ind,:) = flairint_values;
        ventdistfeats(ind,:) = ventdist_values;
        lesvolfeats(ind,:) = lesvol_values;
    end
    
