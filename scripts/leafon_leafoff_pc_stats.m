%% specify input year and paths for las, dtm files

% specify coordinates of study area
aoi = [669640 670000 258910 259250]; % Laegern to fit all polygons
%aoi = [669660 669960 258910 259210]; % Laegern from Felix's script


%% 2010

year = '2010';

% 2010 Spring from lidarlab/data/Laegeren/Spring_2010/Pointcloud/
dir_raw_off = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/Pointcloud/';

% 2010 Summer from lidarlab/data/Laegeren/Summer_2010/Pointcloud/
dir_raw_on = '/Users/scholl/geo_uzh/data/Laegeren/2010_on/Pointcloud/';

% DTM 
dtmFilepath = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/DTM/dtm.mat';


%% 2014 

year = '2014'; 

% 2014 Spring from lidarlab/data/KantonAargau/LeafOff/LAS/
dir_raw_off = '/Users/scholl/geo_uzh/data/Laegeren/2014_off/Pointcloud/';

% 2014 Summer from lidarlab/data/KantonAargau/LeafOn/LAS/
dir_raw_on = '/Users/scholl/geo_uzh/data/Laegeren/2014_on/Pointcloud/';

% DTM
dtmFilepath = '/Users/scholl/geo_uzh/data/Laegeren/2014_off/669000_258000_leaf_off.mat';


%% 2014 Laegeren data in tile form from Kanton Aargau

% search directoy for laz files. if a corresponding las file does not
% exist, convert to laz to als using lasmerge. 
cd(dir_raw_off)
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


cd(dir_raw_on)
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


%% Merge and read las files within study region

raw_off_all = getrawlas(dir_raw_off,aoi);
% raw_off_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2014_off/Pointcloud/data.las');
% raw_off_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2010_off/Pointcloud/data.las');

raw_on_all = getrawlas(dir_raw_on,aoi);
% raw_on_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2014_on/Pointcloud/data.las');
% raw_on_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2010_on/Pointcloud/data.las');


%% interpolate dtm to the raw point coordinates
load(dtmFilepath)
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);

% normalize PC to DTM
raw_off_all.tz = raw_off_all.z;
raw_off_all.z = raw_off_all.tz - interp2(dtm.X,dtm.Y,dtm.z,raw_off_all.x,raw_off_all.y);

raw_on_all.tz = raw_on_all.z;
raw_on_all.z = raw_on_all.tz - interp2(dtm.X,dtm.Y,dtm.z,raw_on_all.x,raw_on_all.y);

% remove noisy values above 50m. 
i = raw_off_all.z<50; 
raw_off = subsetraw(raw_off_all,i);

i = raw_on_all.z<50; 
raw_on = subsetraw(raw_on_all,i);

% plot 
figure; myscatter3(raw_off.x,raw_off.y,raw_off.z,raw_off.z,parula);
title([year ' Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on

figure; myscatter3(raw_on.x,raw_on.y,raw_on.z,raw_on.z,parula);
title([year ' Leaf Off PC Normalized']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on


%% save DMT-normalized, noise-filtered PC mat files

save(['/Users/scholl/geo_uzh/output/laegern/raw_' year '_leafoff_filtered.mat'], 'raw_off');
save(['/Users/scholl/geo_uzh/output/laegern/raw_' year '_leafon_filtered.mat'], 'raw_on'); 


%%  Read in crown polygons
run crown_polygons_laegeren.m


%% create empty array to store stats
n_trees = size(laegernTreeTable_final,1);
% n_trees=20; % for testing, not looping through all trees


% leaf off 
% create structure field for polygon ID, species ID
stats_off.idField = laegernTreeTable_final.idField;
stats_off.species = laegernTreeTable_final.species;

%figure; hold on; title('Leaf Off - Crown Polygons and associated points');

for j = 1:n_trees     
    
    % find raw las points within current polygon
    xpoly = laegernTreeTable_final.xPoly{j};
    ypoly = laegernTreeTable_final.yPoly{j};
    in = inpolygon(raw_off.x,raw_off.y,xpoly,ypoly);
    zpoly = raw_off.z(in);
    
%     % plot current polygon and points within it
%     myscatter3(raw_off.x(in),raw_off.y(in),zpoly,zpoly,parula);
%     z = repmat(45,size(xpoly));
%     d = patch(xpoly,ypoly,z,[0 0.4 0.8]);
%     d.FaceAlpha = 0.4;
    
    % calculate point cloud statistics 
    zpoly_above3m = zpoly(zpoly>3);       % echo heights > 3m above ground
    stats_off.zMax(j,1) = max(zpoly_above3m);       % max height
    stats_off.zMedian(j,1) = median(zpoly_above3m); % median height
    stats_off.zMean(j,1) = mean(zpoly_above3m);     % mean height
    stats_off.zStd(j,1) = std(zpoly_above3m);       % standard deviation
    
    % fraction of single echos 
    r = raw_off.rnnr(in);       
    r1 = ismember(r,11) & (zpoly>3);         
    stats_off.singleEchoFraction(j,1) = sum(r1) / numel(zpoly_above3m);
    
    % fraction of ground echos (single returns < 0.5m)
    g1 = ismember(r,11) & (zpoly<0.5);
    %stats_off.groundEchoFraction(j,1) = sum(g1) / numel(zpoly);
    stats_on.groundEchoFraction(j,1) = sum(g1) / numel(zpoly_above3m);
end


% leaf on
stats_on.idField = laegernTreeTable_final.idField;
stats_on.species = laegernTreeTable_final.species;

for j = 1:n_trees     
     % find raw las points within current polygon
    xpoly = laegernTreeTable_final.xPoly{j};
    ypoly = laegernTreeTable_final.yPoly{j};
    in = inpolygon(raw_on.x,raw_on.y,xpoly,ypoly);
    zpoly = raw_on.z(in);
    
    % calculate point cloud statistics 
    zpoly_above3m = zpoly(zpoly>3);       % echo heights > 3m above ground
    stats_on.zMax(j,1) = max(zpoly_above3m);       % max height
    stats_on.zMedian(j,1) = median(zpoly_above3m); % median height
    stats_on.zMean(j,1) = mean(zpoly_above3m);     % mean height
    stats_on.zStd(j,1) = std(zpoly_above3m);       % standard deviation
    
    % fraction of single echos 
    r = raw_on.rnnr(in);       
    r1 = ismember(r,11) & (zpoly>3);         
    stats_on.singleEchoFraction(j,1) = sum(r1) / numel(zpoly_above3m);
    
    % fraction of ground echos (single returns < 0.5m)
    g1 = ismember(r,11) & (zpoly<0.5);
    %stats_on.groundEchoFraction(j,1) = sum(g1) / numel(zpoly);
    stats_on.groundEchoFraction(j,1) = sum(g1) / numel(zpoly_above3m);
end


%% save stats structs to file

save(['/Users/scholl/geo_uzh/output/laegern/stats_' year '_leafoff.mat'], 'stats_off');
save(['/Users/scholl/geo_uzh/output/laegern/stats_' year '_leafon.mat'], 'stats_on');

%% load lidar mat files 
raw_off_2010 = load('/Users/scholl/geo_uzh/output/laegern/raw_2010_leafoff_filtered.mat');
raw_on_2010 = load('/Users/scholl/geo_uzh/output/laegern/raw_2010_leafon_filtered.mat');
raw_off_2014 = load('/Users/scholl/geo_uzh/output/laegern/raw_2014_leafoff_filtered.mat');
raw_on_2014 = load('/Users/scholl/geo_uzh/output/laegern/raw_2014_leafon_filtered.mat');

raw_2010.leafOff = raw_off_2010.raw_off;
raw_2010.leafOn = raw_on_2010.raw_on;
raw_2014.leafOff = raw_off_2014.raw_off;
raw_2014.leafOn = raw_on_2014.raw_on;

%% restrict to the 8 species in Felix's paper

% load statistic mat files
stats_off_2010 = load('/Users/scholl/geo_uzh/output/laegern/stats_2010_leafoff.mat');
stats_on_2010 = load('/Users/scholl/geo_uzh/output/laegern/stats_2010_leafon.mat'); 
stats_off_2014 = load('/Users/scholl/geo_uzh/output/laegern/stats_2014_leafoff.mat'); 
stats_on_2014 = load('/Users/scholl/geo_uzh/output/laegern/stats_2014_leafon.mat'); 

stats_2010.leafOff = stats_off_2010.stats_2010_leafoff;
stats_2010.leafOn = stats_on_2010.stats_2010_leafon;
stats_2014.leafOff = stats_off_2014.stats_2014_leafoff;
stats_2014.leafOn = stats_on_2014.stats_2014_leafon;



%% keep only species 11 14 22 23 29 31 56 59, each has 40 or more polygons
k = ismember(stats_off_2010.stats_off.species,[11 14 22 23 29 31 56 59]); 
% stats_off_2010 = subsetraw(stats_off_2010.stats_off,k);
% stats_on_2010 = subsetraw(stats_on_2010.stats_on,k);
% stats_off_2014 = subsetraw(stats_off_2014.stats_off,k);
% stats_on_2014 = subsetraw(stats_on_2014.stats_on,k);

stats_off_2010 = subsetraw(stats_2010.leafOff,k);
stats_on_2010 = subsetraw(stats_2010.leafOn,k);
stats_off_2014 = subsetraw(stats_2014.leafOff,k);
stats_on_2014 = subsetraw(stats_2014.leafOn,k);

%% boxplots 


% 2010 
figure('Name','Laegeren 2010 Leaf Off - Leaf On Difference Boxplots'); 

subplot(6,1,1);
dif_zMax = stats_off_2010.zMax - stats_on_2010.zMax; 
boxplot(dif_zMax,stats_off_2010.species);title('max height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,2);
dif_zMedian =  stats_off_2010.zMedian - stats_on_2010.zMedian; 
boxplot(dif_zMedian,stats_off_2010.species);title('median height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-4,3])

subplot(6,1,3);
dif_zStd = stats_off_2010.zStd - stats_on_2010.zStd; 
boxplot(dif_zStd,stats_off_2010.species);title('std height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,4);
dif_zMean = stats_off_2010.zMean - stats_on_2010.zMean; 
boxplot(dif_zMean,stats_off_2010.species);title('mean height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,5);
dif_singleEchoFraction = stats_off_2010.singleEchoFraction - stats_on_2010.singleEchoFraction; 
boxplot(dif_singleEchoFraction,stats_off_2010.species);title('fraction of single echos'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,6);
dif_groundEchoFraction = stats_off_2010.groundEchoFraction - stats_on_2010.groundEchoFraction; 
boxplot(dif_groundEchoFraction,stats_off_2010.species);title('fraction of ground echos'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-0.2,0.4]);


% 2014 
figure('Name','Laegeren 2014 Leaf Off - Leaf On Difference Boxplots'); 

subplot(6,1,1);
dif_zMax = stats_off_2014.zMax - stats_on_2014.zMax; 
boxplot(dif_zMax,stats_off_2014.species);title('max height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,2);
dif_zMedian =  stats_off_2014.zMedian - stats_on_2014.zMedian; 
boxplot(dif_zMedian,stats_off_2014.species);title('median height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,3);
dif_zStd = stats_off_2014.zStd - stats_on_2014.zStd; 
boxplot(dif_zStd,stats_off_2014.species);title('std height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,4);
dif_zMean = stats_off_2014.zMean - stats_on_2014.zMean; 
boxplot(dif_zMean,stats_off_2014.species);title('mean height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,5);
dif_singleEchoFraction = stats_off_2014.singleEchoFraction - stats_on_2014.singleEchoFraction; 
boxplot(dif_singleEchoFraction,stats_off_2014.species);title('fraction of single echos'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,6);
dif_groundEchoFraction = stats_off_2014.groundEchoFraction - stats_on_2014.groundEchoFraction; 
boxplot(dif_groundEchoFraction,stats_off_2014.species);title('fraction of ground echos'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-0.2,0.4]);


%% testing different ground echo fraction calculations 


