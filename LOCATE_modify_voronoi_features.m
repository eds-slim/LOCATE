function [flairfeatmats_modified,ventdistfeatmats_modified,lesvolfeatmats_modified] = LOCATE_modify_voronoi_features(flairfeatmats,ventdistfeatmats,lesvolfeatmats)

flairfeatmats_modified = cell(1,size(flairfeatmats,1));
for nums = 1:size(flairfeatmats,1)
    feats = flairfeatmats{nums};
    feats_modified = [];
    for length = 1:size(feats,1)
        no_of_feats = size(feats,2)/20;
        feature1 = zeros(1,size(feats,2));
        for no_feats = 1:no_of_feats
            feature11 = feats(length,20*(no_feats-1)+1:20*no_feats);
            pos1 = find(feature11 == -1, 1, 'first');
            if isempty(pos1)
                feature1(20*(no_feats-1)+1:20*no_feats) = feature11;
                continue;
            elseif pos1 == 1
                feature11 = 0.*feature11;
            else
                value = feature11(pos1-1);
                feature11(pos1:end) = value;
            end 
            feature1(20*(no_feats-1)+1:20*no_feats) = feature11;
        end            
        feats_modified = [feats_modified;feature1];
    end
    flairfeatmats_modified{nums} = feats_modified;
end


ventdistfeatmats_modified = cell(1,size(flairfeatmats,1));
for nums = 1:size(flairfeatmats,1)
    feats = ventdistfeatmats{nums};
    feats_modified = [];
    for length = 1:size(feats,1)
        feature1 = feats(length,1:size(feats,2));
        pos1 = find(feature1 == -1, 1, 'first');
        if isempty(pos1)
            feats_modified = [feats_modified;feature1];
            continue;
        elseif pos1 == 1
            feature1 = zeros(1,size(feats,2));
        else
            value = feature1(pos1-1);
            feature1(pos1:end) = value;
        end
        feats_modified = [feats_modified;feature1];
    end
    ventdistfeatmats_modified{nums} = feats_modified;
end


lesvolfeatmats_modified = cell(1,size(flairfeatmats,1));
for nums = 1:size(flairfeatmats,1)
    feats = lesvolfeatmats{nums};
    feats_modified = [];
    for length = 1:size(feats,1)
        feature1 = feats(length,1:size(feats,2));
        pos1 = find(feature1 == -1, 1, 'first');
        if isempty(pos1)
            feats_modified = [feats_modified;feature1];
            continue;
        elseif pos1 == 1
            feature1 = zeros(1,size(feats,2));
        else
            value = feature1(pos1-1);
            feature1(pos1:end) = value;
        end 
        feats_modified = [feats_modified;feature1];
    end
    lesvolfeatmats_modified{nums} = feats_modified;
end
end