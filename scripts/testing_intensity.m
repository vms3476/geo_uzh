
base = '/Users/scholl/geo_uzh/';
points = [base 'data/KantonAargau/LeafOff/'];
file = '631000_264000.las';
las_name = [points file];

% need to convert laz to las before using LASread function
if exist(las_name, 'file') == 0
    laszip = '/Users/scholl/LAStools/bin/laszip';
    laz_name = [las_name(1:end-1) 'z'];
    unix([laszip ' -i ' laz_name ' -o ' las_name]);
end

pc = LASread([las_name]);
x = pc.record.x;
y = pc.record.y;
z = pc.record.z;

%%  create a logical indexing variable for all classes except noise points
idxl_noise = ismember(pc.record.classification, 7);

% first return points 
idxl_first = ismember(pc.record.return_number, 1);


%% plot non-noise points colored by return intensity:
figure
scatter3(pc.record.x(~idxl_noise), ...
pc.record.y(~idxl_noise), ...
pc.record.z(~idxl_noise), ...
6, ...
pc.record.intensity(~idxl_noise), ...
'Marker', '.');
colorbar
axis equal tight vis3d
title('Return intensity')
xlabel('X')
ylabel('Y')
ylabel('Z')

%% plot return intensities greater than 400
% plot non-noise points colored by return intensity:

idxl_int = (pc.record.intensity > 200);

figure
scatter3(pc.record.x(~idxl_noise & idxl_int), ...
pc.record.y(~idxl_noise & idxl_int), ...
pc.record.z(~idxl_noise & idxl_int), ...
6, ...
pc.record.intensity(~idxl_noise & idxl_int), ...
'Marker', '.');
colorbar
axis equal tight vis3d
title('Return intensity > 200')
xlabel('X')
ylabel('Y')
ylabel('Z')

%% plot points classified as ground (2)

idxl_grnd = ismember(pc.record.classification, 2);

figure
scatter3(pc.record.x(idxl_grnd & idxl_first), ...
pc.record.y(idxl_grnd & idxl_first), ...
pc.record.z(idxl_grnd & idxl_first), ...
6, ...
pc.record.classification(idxl_grnd & idxl_first), ...
'Marker', '.');
colorbar
axis equal tight vis3d
title('Classified as Ground')
xlabel('X')
ylabel('Y')
ylabel('Z')

%% plot first returns classified as vegetation (3,4,5) 

class = pc.record.classification;
idxl_veg = ismember(class,3) | ismember(class,4) | ismember(class,5);

figure
hold on;
scatter3(pc.record.x(idxl_veg & idxl_first), ...
pc.record.y(idxl_veg & idxl_first), ...
pc.record.z(idxl_veg & idxl_first), ...
6, ...
pc.record.z(idxl_veg & idxl_first), ...
'Marker', '.');
colorbar
axis equal tight vis3d
title('First Returns Classified as Veg')
xlabel('X')
ylabel('Y')
ylabel('Z')

% now plot the first returns not vegetation
scatter3(pc.record.x(~idxl_veg & idxl_first), ...
pc.record.y(~idxl_veg & idxl_first), ...
pc.record.z(~idxl_veg & idxl_first), ...
'Marker', '.');

%% create shp file to convert to kml and be displayed on map.geo.admin.ch

outputFilepath = [las_name(1:end-4) '.shp'];
extent = LASextent([points '631000_264000/'], outputFilepath,'fig', true, 'verbose', true);
shp = shaperead(outputFilepath);

%%

[models, refmat] = elevationModels([x, y, z], ...
  class, ...
  'classTerrain', [2], ...
  'cellResolution', 1, ...
  'maxFillArea', inf, ...
  'smoothingFilter', [], ...
  'outputModels', {'terrain', 'surface', 'height', 'density'}, ...
  'fig', true, ...
  'verbose', true);

  spatialRef = refmatToMapRasterReference(refmat, size(models.height.values), ...
    'rasterInterpretation', 'cells');
