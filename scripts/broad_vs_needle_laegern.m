% testing deciduous vs. coniferous identification at Laegarn site 

load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat');

% isolate the column containing species data 
species = laegernTreeTable_final.species;

% specify the code numbers that refer to coniferous vs. deciduousleaf trees
code_coniferous = [11, 14, 15];
code_deciduous = [21, 22, 23, 28, 29, 31, 41, 51, 56, 59];

% find indices within tree table of coniferous vs. deciduousleaf trees
idx_coniferous = sum(species(:)== code_coniferous,2);
data_idx_coniferous = find(idx_coniferous==1);

idx_deciduous = sum(species(:)== code_deciduous,2);
data_idx_deciduous = find(idx_deciduous==1);


% Plot deciduousleaf tree polygons
figure; 
for i = 1:numel(data_idx_deciduous)
    x = laegernTreeTable_final.xPoly{data_idx_deciduous(i)};
    y = laegernTreeTable_final.yPoly{data_idx_deciduous(i)};
    d = patch(x,y,[0 0.4 0.8]);
end

for i = 1:numel(data_idx_coniferous)
    x = laegernTreeTable_final.xPoly{data_idx_coniferous(i)};
    y = laegernTreeTable_final.yPoly{data_idx_coniferous(i)};
    c = patch(x,y,[0, 0.8, 0.4]);
end

legend([d,c],'Deciduous','Coniferous')
xt = get(gca, 'XTick');
set(gca, 'FontSize', 14)
xlabel('Easting [m]','FontSize',16); xtickangle(45)
ylabel('Northing [m]','FontSize',16); ytickangle(45)
title('Largern Crown Map - Deciduous / Coniferous','FontSize',16)

%% Read Laegeren leaf off data and isolate deciduous trees 




