
% define input files and parameters
base = '/Users/scholl/geo_uzh/';
points = [base 'data/KantonAargau/LeafOff/'];
file = '631000_264000.las';
las_name = [points file];
res = 0.5;
cmap = parula; 

% need to convert laz to las before using LASread function
if exist(las_name, 'file') == 0
    laszip = '/Users/scholl/LAStools/bin/laszip';
    laz_name = [las_name(1:end-1) 'z'];
    unix([laszip ' -i ' laz_name ' -o ' las_name]);
end

% read in point cloud data
pc = LASread(las_name);
x = pc.record.x;
y = pc.record.y;
z = pc.record.z;
classification = pc.record.classification;


% clip subsets for deciduous and coniferous areas

% deciduous
xv = [631760; 631760; 631860; 631860];
yv = [264610; 264710; 264710; 264610];
outputFilepath = [las_name(1:end-4) '_subset.las'];
pc_s_deciduous = LASclip(pc, [xv, yv], outputFilepath, 'verbose', true);

% deciduous subset
xv = [631760; 631760; 631860; 631860];
yv = [264610; 264710; 264710; 264610];
outputFilepath = [las_name(1:end-4) '_subset_deciduous.las'];
s = LASclip(pc, [xv, yv], outputFilepath, 'verbose', true);
s_deciduous = LASread(outputFilepath);

% plot subset
figure; 
scatter3(s_deciduous.record.x(:), s_deciduous.record.y(:), ...
         s_deciduous.record.z(:), 6, s_deciduous.record.z(:), ...
         'Marker', '.');
colorbar; axis equal tight vis3d
title('Deciduous Subset'); xlabel('X'); ylabel('Y'); ylabel('Z')

% coniferous subset
xv = [631840; 631840; 631940; 631940];
yv = [264300; 264400; 264400; 264300];
outputFilepath = [las_name(1:end-4) '_subset_coniferous.las'];
s = LASclip(pc, [xv, yv], outputFilepath, 'verbose', true);
s_coniferous = LASread(outputFilepath);

% plot subset
figure; 
scatter3(s_coniferous.record.x(:), s_coniferous.record.y(:), ...
         s_coniferous.record.z(:), 6, s_coniferous.record.z(:), ...
         'Marker', '.');
colorbar; axis equal tight vis3d
title('Coniferous Subset'); xlabel('X'); ylabel('Y'); ylabel('Z')

%% dtm
[models, refmat] = elevationModels([x, y, z], ...
  classification, ...
  'classTerrain', [2], ...
  'cellResolution', res, ...
  'maxFillArea', inf, ...
  'smoothingFilter', [], ...
  'outputModels', {'terrain', 'surface', 'height', 'density'}, ...
  'fig', true, ...
  'verbose', true);

spatialRef = refmatToMapRasterReference(refmat, size(models.height.values), ...
            'rasterInterpretation', 'cells');
        
        
geotiffwrite([outputFilepath(1:end-4) '.tif'], ...
            single(models.height.values), spatialRef, ...
            'CoordRefSysCode', 'EPSG:21781');

%create grid of Easting / Northing coordinates to plot dtm

grid.x = unique(models.terrain.interpolant.Points(:,1)');

grid.y = unique(models.terrain.interpolant.Points(:,2)');
[grid.X, grid.Y]  = meshgrid(grid.x,grid.y);
figure; 
scatter3(grid.X(:), grid.Y(:), ...
         models.terrain.values(:), 6, models.terrain.values(:), ...
         'Marker', '.');
colorbar; axis equal tight vis3d
title('DTM'); xlabel('X'); ylabel('Y'); ylabel('Z')

%% read dtm from geotiff
dtm = geotiffread


%%    

% index based on first returns only
i = ismember(pc.record.classification, 2);

x1 = x(i);
y1 = y(i);
z1 = z(i);

%% create dtm raster 

dxy = res;
xv = min(x1):dxy:max(x1);
yv = min(y1):dxy:max(y1);
zv = min(z1):dxy:max(z1);
xv = [xv, xv(end)+ceil(mod(max(x), dxy))*dxy];
yv = [yv, yv(end)+ceil(mod(max(y), dxy))*dxy];
zv = [zv, zv(end)+ceil(mod(max(z), dxy))*dxy];
[~, sub, raster] = rasterize([x1, y1, z1], xv, yv, zv, z1, @(x) numel(x));

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


%% isolate first returns

%subtract first return DSM - DTM
















%% plot mparkan dtm

figure
gh11 = imagesc(models.terrain.values);
colormap(gray)
set(gh11, 'AlphaData', models.mask)
title('Terrain Model')
axis equal tight
colorbar


[grid.X,grid.Y] = meshgrid(unique(models.terrain.interpolant.Points(:,1)),...
                           unique(models.terrain.interpolant.Points(:,2)));
                       
dtm_interp = interp2(grid.X, grid.Y,...
                     models.terrain.values,...
                     grid.x(:), grid.y(:));
figure;
xmax = xmin + 1000;
ymax = ymin + 1000;
X = xmin:res:xmax;
Y = ymin:res:ymax;
[dtmX, dtmY] = meshgrid(X,Y);
myscatter3(dtmX(:), dtmY(:), models.terrain.values(:), models.terrain.values(:), cmap);
myscatter3(models.terrain.interpolant.Points(:,1)',...
           models.terrain.interpolant.Points(:,2)',...
           models.terrain.interpolant.Values(:)',...
           models.terrain.interpolant.Values(:)', cmap);
title('DTM'); swisstick




