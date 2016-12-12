function [  ] = overlay_polygons( laegernTreeTable_final )

% load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');

species = laegernTreeTable_final.species;

k = ismember(species,[11 14 22 23 29 31 56 59]); 
trees = laegernTreeTable_final(k,:);

% specify the code numbers that refer to coniferous vs. deciduousleaf trees
code_coniferous = [11, 14];
code_deciduous = [22, 23, 29, 31, 56, 59];

% find indices within tree table of coniferous vs. deciduousleaf trees
idx_coniferous = sum(species(:)== code_coniferous,2);
data_idx_coniferous = find(idx_coniferous==1);

idx_deciduous = sum(species(:)== code_deciduous,2);
data_idx_deciduous = find(idx_deciduous==1);


for i = 1:numel(data_idx_deciduous)
    x = laegernTreeTable_final.xPoly{data_idx_deciduous(i)};
    y = laegernTreeTable_final.yPoly{data_idx_deciduous(i)};
    z = repmat(50,size(x));
    d = patch(x,y,z,[0 0.4 0.8],'FaceAlpha',0,'EdgeColor',[0.2 0.2 1]);
end

for i = 1:numel(data_idx_coniferous)
    x = laegernTreeTable_final.xPoly{data_idx_coniferous(i)};
    y = laegernTreeTable_final.yPoly{data_idx_coniferous(i)};
    z = repmat(50,size(x));
    c = patch(x,y,z,[0, 0.8, 0.4],'FaceAlpha',0,'EdgeColor',[0 0.8 0.4]);
end


end

