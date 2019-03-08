function [lesionmask, index_mask, index_numbers] =  LOCATE_Voronoi_tessellation_resized(lesionmask, biancamask,inv_factor)
% Function for extracting test features for LOCATE
%   Copyright - FMRIB, WIN, University of Oxford
%   Vaanathi Sundaresan - 25/05/2018
%
%   lesionmask - BIANCA lesion probability maps
%   biancamask - WM matter with SC structures used in BIANCA
%   Inv_factor - inverse of up/downsampling factor to resample images to
%   original dimensions

    % Smoothing the images with Gaussian kernel to remove spurious lesion
    % voxels
    lesionmask_smooth = imgaussfilt3(lesionmask,1,'FilterSize',5);
    
    % Determine local maxima of the lesion mask
    regmax = imregionalmax(lesionmask_smooth);
    regmax = regmax.*(biancamask>0);
    regmax_props = regionprops(bwlabeln(regmax),'Centroid');
    temp = struct2cell(regmax_props);
    regmax_centroid = cell2mat(temp(1,:)');
    
    % Create dummy points to restrict the voronoi points from extending to
    % infinity
    regmax_centroid_dummy = [1,1,1;size(lesionmask,1),1,1;1,size(lesionmask,2),1;size(lesionmask,1),size(lesionmask,2),1;...
    1,1,size(lesionmask,3);size(lesionmask,1),1,size(lesionmask,3);1,size(lesionmask,2),size(lesionmask,3);...
    size(lesionmask,1),size(lesionmask,2),size(lesionmask,3)];
    regmax_centroid= [regmax_centroid_dummy;regmax_centroid];
    
    % Voronoi Tessellation
    [vertices,indices] = voronoin(regmax_centroid);
    index_mask = zeros(size(lesionmask));
    voronoi_area = zeros(1,size(indices,1)-size(regmax_centroid_dummy,1));
    
    % Iterating over voronoi regions
    for ind = size(regmax_centroid_dummy,1)+1:size(indices,1) %        
        temp = zeros(size(lesionmask));
        reqd_indices = indices{ind};
        reqd_indices(reqd_indices == 1) = [];
        reqd_vertices = floor(vertices(reqd_indices,:));
        reqd_vertices1 = reqd_vertices(:,1);
        reqd_vertices2 = reqd_vertices(:,2);
        reqd_vertices3 = reqd_vertices(:,3);

        reqd_vertices1 = max(reqd_vertices1,ones(numel(reqd_vertices1),1));
        reqd_vertices1 = min(reqd_vertices1,size(lesionmask,2).*ones(numel(reqd_vertices1),1));
        reqd_vertices2 = max(reqd_vertices2,ones(numel(reqd_vertices1),1));
        reqd_vertices2 = min(reqd_vertices2,size(lesionmask,1).*ones(numel(reqd_vertices1),1));
        reqd_vertices3 = max(reqd_vertices3,ones(numel(reqd_vertices1),1));
        reqd_vertices3 = min(reqd_vertices3,size(lesionmask,3).*ones(numel(reqd_vertices1),1));
        reqd_vertices = [reqd_vertices1,reqd_vertices2,reqd_vertices3];
        % Getting voronoi vertices
        for num = 1:numel(reqd_indices)
            temp(reqd_vertices(num,2),reqd_vertices(num,1), reqd_vertices(num,3)) = 1;
        end
        
        % And forming voronoi regions (as convex hulls)
        voronoi_region_mask = bwconvhull3d(temp);
        voronoi_region_mask = voronoi_region_mask & biancamask;
        voronoi_region_mask(index_mask>0) = 0;
        index = voronoi_region_mask.*(ind-size(regmax_centroid_dummy,1));
        % Forming Voronoi index maps
        index_mask = index_mask + index;
        voronoi_area(ind-size(regmax_centroid_dummy,1)) = sum(voronoi_region_mask(:));
    end
%     indexmask = index_mask;
%     % Resizing the images back to the original dimensions
%     lesionmask = imresizen(lesionmask,inv_factor);
%     biancamask = imresizen(single(biancamask),inv_factor);
%     index_mask = imresizen(index_mask,inv_factor,'nearest');
%     index_mask = round(index_mask);
%     index_mask = index_mask.*(biancamask);
    index_numbers = setdiff(union(index_mask(:),[]),0);
    
end

       