
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

% create dtm




% isolate first returns

% subtract first return DSM - DTM







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




