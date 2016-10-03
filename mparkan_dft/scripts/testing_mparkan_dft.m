% testing_mparkan_dft code scripts using example syntax within the provided
% documentation of the Digital Forestry Toolbox by M. markan

points = '/media/Data/victoria/geo_uzh/mparkan_dft/data/measurements/vector/als/';
pc = LASread([points 'zh_2014_coniferous.las']);
x = pc.record.x;
y = pc.record.y;
z = pc.record.z;

% create a logical indexing variable for all classes except noise points
idxl_noise = ismember(pc.record.classification, 7);

% plot non-noise points colored by return intensity:
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

%% export extent of LAS files to a single ESRI shapefile

% there was an issue with the code adding a backslash onto the end of the
% file path! 
outputFilepath = [points 'zh_2014_coniferous.shp'];
extent = LASextent(points, outputFilepath,'fig', true, 'verbose', true);

%% Extract a 2D cross-section from a 3D point cloud

 width = 2;
 p0 = [699501, 271206];
 p1 = [699698, 271199];
 [index, footprint, profile] = crossSection([x, y, z], width, p0, p1, 'verbose', false, 'fig', true);

%% Convert a 3D point cloud to a 2D/3D raster

 dxy = 0.5;
 xv = min(x):dxy:max(x);
 yv = min(y):dxy:max(y);
 zv = min(z):dxy:max(z);
 xv = [xv, xv(end)+ceil(mod(max(x), dxy))*dxy];
 yv = [yv, yv(end)+ceil(mod(max(y), dxy))*dxy];
 zv = [zv, zv(end)+ceil(mod(max(z), dxy))*dxy];
 [~, sub, raster] = rasterize([x, y, z], xv, yv, zv, z, @(x) numel(x));

 nrows = length(yv);
 ncols = length(xv);
 nlevels = length(zv);

 ind = sub2ind([nrows, ncols, nlevels], sub(:,2), sub(:,1), sub(:,3)); % 3d (row, col, lev) linear index for each point
 [idxn_cells, ~, idxn_voxel_to_point] = unique(ind); % assign raster cell metrics to corresponding points
 [row_voxel, col_voxel, lev_voxel] = ind2sub([nrows, ncols, nlevels], idxn_cells);
 density = raster(idxn_cells);

 figure
 scatter3(col_voxel, row_voxel, lev_voxel, 20, ...
     density, ...
     'Marker', '.');
 xlabel('col')
 ylabel('row')
 zlabel('level')
 title('Point count per voxel')
 colorbar
 axis equal tight vis3d

 %% Compute elevation models (i.e. terrain, surface, height) from a 
 %  3D classified point cloud
 
 classification = pc.record.classification;

 [models, refmat] = elevationModels([x, y, z], ...
  classification, ...
  'classTerrain', [2], ...
  'classSurface', [18], ...
  'cellResolution', 1, ...
  'maxFillArea', inf, ...
  'smoothingFilter', [], ...
  'outputModels', {'terrain', 'surface', 'height', 'density'}, ...
  'fig', true, ...
  'verbose', true);

  spatialRef = refmatToMapRasterReference(refmat, size(models.height.values), ...
    'rasterInterpretation', 'cells');

% original command doesn't have the first quote for the output path
  geotiffwrite('/media/Data/victoria/geo_uzh/mparkan_dft/data/measurements/raster/als/so_2014_woodland_pasture.tif', ...
    single(models.height.values), spatialRef, ...
    'CoordRefSysCode', 'EPSG:21781');

%% Computes the echo ratio of a 3D point cloud

 echo_ratio = laserEchoRatio([x, y, z], 'rasterResolution', 1, 'verbose', true, 'fig', true);

     
%% Determines the pulse number associated with the individual returns 
%  (sorted by acquisition GPS time)

return_number = 1;
[pulse_number, pulse_returns, ind_pulse] = laserPulses(return_number);

%% Tree metrics %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract individual tree crowns from a 3D point cloud using a modified 
% version of the top down region growing algorithm described in Li et al. (2012)

xyz = [pc.record.x, pc.record.y, pc.record.z];
classification = pc.record.classification;

[point, tree] = treeRegionGrowing([x y z], ...
classification, ...
'classTerrain', [2], ...
'classHighVegetation', [4,5,12], ...
'coordinateResolution', 0.5, ...
'normalizeElevation', true, ...
'peakSearchRadius', 1.5, ...
'minPeakSpacing', [0 1.5; 15 3], ...
'minSampleSize', 20, ...
'minSamplingRadius', 5, ...
'maxSamplingRadius', 18, ...
'verbose', true, ...
'fig', true);


%% Extract individual tree crowns from a raster Canopy Height Model (CHM) 

info = geotiffinfo('/media/Data/victoria/geo_uzh/mparkan_dft/data/measurements/raster/chm/so_2014_woodland_pasture.tif');
[chm, ~] = geotiffread('/media/Data/victoria/geo_uzh/mparkan_dft/data/measurements/raster/chm/so_2014_woodland_pasture.tif');


% Find local maxima (peaks) coordinates in a raster Canopy Height Model (CHM)

[crh, xyh] = treePeaks(chm, info.RefMatrix, ...
 'method', 'fixedRadius', ...
 'windowRadius', 3, ...
 'fig', true);

% Marker-controlled watershed segmentation described in Kwak et al. (2007)

labels = treeWatershed(chm, crh, ...
 'minHeight', 1, ...
 'fig', true);

% Computes a canopy crown cover map using either a convolution window or 
% the Delaunay triangulation method described in Eysn et al. (2012)

cover = canopyCover(chm, info.RefMatrix, ...
    'method', 'triangulation', ...
    'stems', crh, ...
    'crowns', labels, ...
    'verbose', true, ...
    'fig', true);
 
