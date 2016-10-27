%% specify input paths for las, dtm files

aoi = [669640 670000 258910 259250]; % Laegern to fit all polygons

%%  read laegern las files

disp('Reading las files...');

% 2010
dir_pc_2010_off = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/Pointcloud/';
dir_pc_2010_on = '/Users/scholl/geo_uzh/data/Laegeren/2010_on/Pointcloud/';
dir_dtm_2010 = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/DTM/dtm.mat';

% 2014 
dir_pc_2014_off = '/Users/scholl/geo_uzh/data/Laegeren/2014_off/Pointcloud/';
dir_pc_2014_on = '/Users/scholl/geo_uzh/data/Laegeren/2014_on/Pointcloud/';
dir_dtm_2014 = '/Users/scholl/geo_uzh/data/Laegeren/2014_off/669000_258000_leaf_off.mat';


% read las point cloud data

% 2010
pc_2010.leafoff = getrawlas(dir_pc_2010_off,aoi);
pc_2010.leafon = getrawlas(dir_pc_2010_on,aoi);

% 2014
% since in tile form, convert from laz to las before merging
cd(dir_pc_2014_off)
files = dir('*.laz');
if numel(files) > 0 
    for i = 1:numel(files)
        laz_name = [files(i).name];
        las_name = [laz_name(1:end-1) 's'];
        if exist(las_name, 'file') == 0
            laszip = '/Users/scholl/LAStools/bin/laszip';
            unix([laszip ' -i ' laz_name ' -o ' las_name]);
        end
    end
end
cd(dir_pc_2014_on)
files = dir('*.laz');
if numel(files) > 0
    for i = 1:numel(files)
        laz_name = [files(i).name];
        las_name = [laz_name(1:end-1) 's'];
        if exist(las_name, 'file') == 0
            laszip = '/Users/scholl/LAStools/bin/laszip';
            unix([laszip ' -i ' laz_name ' -o ' las_name]);
        end
    end
end
pc_2014.leafoff = getrawlas(dir_pc_2014_off,aoi);
pc_2014.leafon = getrawlas(dir_pc_2014_on,aoi);

%% load dtm   

disp('Normalizing point cloud based on DTM...');

load(dir_dtm_2010); dtm_2010 = dtm; % 2010
load(dir_dtm_2014); dtm_2014 = dtm; % 2014


% interpolate dtm to the raw point coordinates, subtract PC - DTM
    
%2010
[dtm_2010.X,dtm_2010.Y] = meshgrid(dtm_2010.x,dtm_2010.y); 
pc_2010.leafoff.tz = pc_2010.leafoff.z;
pc_2010.leafoff.z = pc_2010.leafoff.tz - ...
    interp2(dtm_2010.X,dtm_2010.Y,dtm_2010.z,pc_2010.leafoff.x,pc_2010.leafoff.y);

pc_2010.leafon.tz = pc_2010.leafon.z;
pc_2010.leafon.z = pc_2010.leafon.tz - ...
    interp2(dtm_2010.X,dtm_2010.Y,dtm_2010.z,pc_2010.leafon.x,pc_2010.leafon.y);

%2014
[dtm_2014.X,dtm_2014.Y] = meshgrid(dtm_2014.x,dtm_2014.y); 
pc_2014.leafoff.tz = pc_2014.leafoff.z;
pc_2014.leafoff.z = pc_2014.leafoff.tz - ...
    interp2(dtm_2014.X,dtm_2014.Y,dtm_2014.z,pc_2014.leafoff.x,pc_2014.leafoff.y);

pc_2014.leafon.tz = pc_2014.leafon.z;
pc_2014.leafon.z = pc_2014.leafon.tz - ...
    interp2(dtm_2014.X,dtm_2014.Y,dtm_2014.z,pc_2014.leafon.x,pc_2014.leafon.y);
    
    
% remove noisy points taller then 50m
pc_2010.leafoff = subsetraw(pc_2010.leafoff,(pc_2010.leafoff.z<50));
pc_2010.leafon = subsetraw(pc_2010.leafon,(pc_2010.leafon.z<50));
pc_2014.leafoff = subsetraw(pc_2014.leafoff,(pc_2014.leafoff.z<50));
pc_2014.leafon = subsetraw(pc_2014.leafon,(pc_2014.leafon.z<50));

    
% plot normalized point cloud

% 2010
figure; myscatter3(pc_2010.leafoff.x,pc_2010.leafoff.y,pc_2010.leafoff.z,pc_2010.leafoff.z,parula);
title(['2010 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on

figure; myscatter3(pc_2010.leafon.x,pc_2010.leafon.y,pc_2010.leafon.z,pc_2010.leafon.z,parula);
title(['2010 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on


%2014
figure; myscatter3(pc_2014.leafoff.x,pc_2014.leafoff.y,pc_2014.leafoff.z,pc_2014.leafoff.z,parula);
title(['2014 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on

figure; myscatter3(pc_2014.leafon.x,pc_2014.leafon.y,pc_2014.leafon.z,pc_2014.leafon.z,parula);
title(['2014 Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on
    

% save DTM-normalized, 50m noise-filtered PC as mat files
pc_2010_leafoff = pc_2010.leafoff; pc_2010_leafon = pc_2010.leafon;
pc_2014_leafoff = pc_2014.leafoff; pc_2014_leafon = pc_2014.leafon;
save('/Users/scholl/geo_uzh/output/laegern/pc_2010_leafoff_filtered.mat', 'pc_2010_leafoff');
save('/Users/scholl/geo_uzh/output/laegern/pc_2010_leafon_filtered.mat', 'pc_2010_leafon');
save('/Users/scholl/geo_uzh/output/laegern/pc_2014_leafoff_filtered.mat', 'pc_2014_leafoff');
save('/Users/scholl/geo_uzh/output/laegern/pc_2014_leafon_filtered.mat', 'pc_2014_leafon');

% %% load DTM-normalized point cloud data
% 
%load('/Users/scholl/geo_uzh/output/laegern/pc_2010_leafoff_filtered.mat')
%load('/Users/scholl/geo_uzh/output/laegern/pc_2010_leafon_filtered.mat')
%load('/Users/scholl/geo_uzh/output/laegern/pc_2014_leafoff_filtered.mat')
%load('/Users/scholl/geo_uzh/output/laegern/pc_2014_leafon_filtered.mat')

%% Read in crown polygons

disp('Reading crown polygons...');

run crown_polygons_laegeren.m

%% Statistics

disp('Computing point cloud statistics...');

% create structure to hold statistics data
stats.idField = laegernTreeTable_final.idField;
stats.species = laegernTreeTable_final.species;
stats.xPoly = laegernTreeTable_final.xPoly;
stats.yPoly = laegernTreeTable_final.yPoly;

n_trees = size(laegernTreeTable_final,1);

for j = 1:n_trees     
    
    j
    
    % find raw las points within current polygon
    xpoly = laegernTreeTable_final.xPoly{j};
    ypoly = laegernTreeTable_final.yPoly{j};
    
    in_2010_off = inpolygon(pc_2010_leafoff.x,pc_2010_leafoff.y,xpoly,ypoly);
    in_2010_on = inpolygon(pc_2010_leafon.x,pc_2010_leafon.y,xpoly,ypoly); 
    in_2014_off = inpolygon(pc_2014_leafoff.x,pc_2014_leafoff.y,xpoly,ypoly); 
    in_2014_on = inpolygon(pc_2014_leafon.x,pc_2014_leafon.y,xpoly,ypoly);
    
    zpoly_2010_off = pc_2010_leafoff.z(in_2010_off);
    zpoly_2010_on = pc_2010_leafon.z(in_2010_on);
    zpoly_2014_off = pc_2014_leafoff.z(in_2014_off);
    zpoly_2014_on = pc_2014_leafon.z(in_2014_on); 
    
    % echo heights > 3m above ground
    zpoly_2010_off_above3m = zpoly_2010_off(zpoly_2010_off>3);
    zpoly_2010_on_above3m = zpoly_2010_on(zpoly_2010_on>3);
    zpoly_2014_off_above3m = zpoly_2014_off(zpoly_2014_off>3);
    zpoly_2014_on_above3m = zpoly_2014_on(zpoly_2014_on>3);
    
    % maximum height
    stats_2010.off.zMax(j,1) = max(zpoly_2010_off_above3m);
    stats_2010.on.zMax(j,1) = max(zpoly_2010_on_above3m);
    stats_2014.off.zMax(j,1) = max(zpoly_2014_off_above3m);
    stats_2014.on.zMax(j,1) = max(zpoly_2014_on_above3m);    
    
    % median height
    stats_2010.off.zMedian(j,1) = median(zpoly_2010_off_above3m);
    stats_2010.on.zMedian(j,1) = median(zpoly_2010_on_above3m);
    stats_2014.off.zMedian(j,1) = median(zpoly_2014_off_above3m);
    stats_2014.on.zMedian(j,1) = median(zpoly_2014_on_above3m);
    
    % mean height  
    stats_2010.off.zMean(j,1) = mean(zpoly_2010_off_above3m);
    stats_2010.on.zMean(j,1) = mean(zpoly_2010_on_above3m);    
    stats_2014.off.zMean(j,1) = mean(zpoly_2014_off_above3m);
    stats_2014.on.zMean(j,1) = mean(zpoly_2014_on_above3m);     
    
    % standard deviation   
    stats_2010.off.zStd(j,1) = std(zpoly_2010_off_above3m);
    stats_2010.on.zStd(j,1) = std(zpoly_2010_on_above3m);
    stats_2014.off.zStd(j,1) = std(zpoly_2014_off_above3m);
    stats_2014.on.zStd(j,1) = std(zpoly_2014_on_above3m);
    
%     % fraction of single echos 
%     r_2010_off = raw_off.rnnr(in);  
%     
%     r1 = ismember(r,11) & (zpoly>3);         
%     stats_off.singleEchoFraction(j,1) = sum(r1) / numel(zpoly_above3m);
%     
%     % fraction of ground echos (single returns < 0.5m)
%     g1 = ismember(r,11) & (zpoly<0.5);
%     %stats_off.groundEchoFraction(j,1) = sum(g1) / numel(zpoly);
%     stats_on.groundEchoFraction(j,1) = sum(g1) / numel(zpoly_above3m);
end

%% save stats structs to file

disp('Saving statistics data to file...');

stats_2010_leafoff = stats_2010.off; stats_2010_leafon = stats_2010.on;
stats_2014_leafoff = stats_2014.off; stats_2014_leafon = stats_2014.on;

save('/Users/scholl/geo_uzh/output/laegern/stats_2010_leafoff.mat', 'stats_2010_leafoff');
save('/Users/scholl/geo_uzh/output/laegern/stats_2010_leafon.mat', 'stats_2010_leafon');
save('/Users/scholl/geo_uzh/output/laegern/stats_2014_leafoff.mat', 'stats_2014_leafoff');
save('/Users/scholl/geo_uzh/output/laegern/stats_2014_leafon.mat', 'stats_2014_leafon');

save('/Users/scholl/geo_uzh/output/laegern/stats.mat', 'stats');


%% Species - keep only 11 14 22 23 29 31 56 59, each has 40 or more polygons
k = ismember(stats.species,[11 14 22 23 29 31 56 59]); 
stats_2010_leafoff = subsetraw(stats_2010_leafoff,k);
stats_2010_leafon = subsetraw(stats_2010_leafon,k);
stats_2014_leafoff = subsetraw(stats_2014_leafoff,k);
stats_2014_leafon = subsetraw(stats_2014_leafon,k);
species = stats.species(k);

%% boxplots

disp('Boxplots...');

% 2010 

figure('Name','Laegeren 2010 Leaf Off - Leaf On Difference Boxplots'); 

subplot(6,1,1);
dif_zMax = stats_2010_leafoff.zMax - stats_2010_leafon.zMax; 
boxplot(dif_zMax,species);title('max height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,2);
dif_zMedian = stats_2010_leafoff.zMedian - stats_2010_leafon.zMedian; 
boxplot(dif_zMedian,species);title('median height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,3);
dif_zStd = stats_2010_leafoff.zStd - stats_2010_leafon.zStd; 
boxplot(dif_zStd,species);title('std height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,4);
dif_zMean = stats_2010_leafoff.zMean - stats_2010_leafon.zMean; 
boxplot(dif_zMean,species);title('mean height'); 
line([0 20],[0 0],'color','k','linewidth',1)

% subplot(6,1,5);
% dif_singleEchoFraction = stats_off_2010.singleEchoFraction - stats_on_2010.singleEchoFraction; 
% boxplot(dif_singleEchoFraction,stats_off_2010.species);title('fraction of single echos'); 
% line([0 20],[0 0],'color','k','linewidth',1)
% 
% subplot(6,1,6);
% dif_groundEchoFraction = stats_off_2010.groundEchoFraction - stats_on_2010.groundEchoFraction; 
% boxplot(dif_groundEchoFraction,stats_off_2010.species);title('fraction of ground echos'); 
% line([0 20],[0 0],'color','k','linewidth',1); ylim([-0.2,0.4]);


% 2014 
figure('Name','Laegeren 2014 Leaf Off - Leaf On Difference Boxplots'); 

subplot(6,1,1);
dif_zMax = stats_2014_leafoff.zMax - stats_2014_leafon.zMax; 
boxplot(dif_zMax,species);title('max height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,2);
dif_zMedian = stats_2014_leafoff.zMedian - stats_2014_leafon.zMedian; 
boxplot(dif_zMedian,species);title('median height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,3);
dif_zStd = stats_2014_leafoff.zStd - stats_2014_leafon.zStd; 
boxplot(dif_zStd,species);title('std height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,4);
dif_zMean = stats_2014_leafoff.zMean - stats_2014_leafon.zMean; 
boxplot(dif_zMean,species);title('mean height'); 
line([0 20],[0 0],'color','k','linewidth',1)
