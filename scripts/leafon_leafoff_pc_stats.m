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
        las_name = laz_name(1:end-1);
        if exist(las_name, 'file') == 0
            laszip = '/Users/scholl/LAStools/bin/laszip';
            unix([laszip ' -i ' laz_name ' -o ' las_name 's']);
        end
    end
end

cd(dir_raw_on)
files = dir('*.laz');
if numel(files) > 0
    for i = 1:numel(files)
        laz_name = [files(i).name];
        las_name = laz_name(1:end-1);
        if exist(las_name, 'file') == 0
            laszip = '/Users/scholl/LAStools/bin/laszip';
            unix([laszip ' -i ' laz_name ' -o ' las_name 's']);
        end
    end
end


%% Merge and read las files within study area

raw_off_all = getrawlas(dir_raw_off,aoi);
%raw_off_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2010_off/Pointcloud/data.las');

raw_on_all = getrawlas(dir_raw_on,aoi);
%raw_on_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2010_on/Pointcloud/data.las');


%% interpolate dtm to the raw point coordinates
load(dtmFilepath)
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);

raw_off_all.tz = raw_off_all.z;
raw_off_all.z = raw_off_all.tz - interp2(dtm.X,dtm.Y,dtm.z,raw_off_all.x,raw_off_all.y);

raw_on_all.tz = raw_on_all.z;
raw_on_all.z = raw_on_all.tz - interp2(dtm.X,dtm.Y,dtm.z,raw_on_all.x,raw_on_all.y);

% remove noisy values above 50m. 
i = raw_off_all.z<50; 
raw_off = subsetraw(raw_off_all,i);

i = raw_on_all.z<50; 
raw_on = subsetraw(raw_on_all,i);

% remove points below a height of 3m, keep only crown returns. 
i = raw_off.z>3; 
raw_off = subsetraw(raw_off,i);

i = raw_on.z>3; 
raw_on = subsetraw(raw_on,i);

%plot 
figure; myscatter3(raw_off.x,raw_off.y,raw_off.z,raw_off.z,parula);
title('Leaf Off PC Normalized'); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on

figure; myscatter3(raw_on.x,raw_on.y,raw_on.z,raw_on.z,parula);
title('Leaf On PC Normalized'); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on


%% PC Statistics for leaf on - leaf off

% subset based on crown polygons
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
    stats_off.zMax(j,1) = max(zpoly);         % max height
    stats_off.zMedian(j,1) = median(zpoly);   % median height
    stats_off.zMean(j,1) = mean(zpoly);       % mean height
    stats_off.zStd(j,1) = std(zpoly);         % standard deviation
    
%     % fraction of single echos = # single / # total echos within crown
%     % fraction of ground echos = # ground / # total echos within crown
%     % compute statistics for the first returns? 
%     r = raw2010_off.rnnr(in);  
%     r1 = sum(r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71);
%     stats.firstEchoFraction(j,1) = r1 / numel(r);
    
end


% leaf on
stats_on.idField = laegernTreeTable_final.idField;
stats_on.species = laegernTreeTable_final.species;

% figure; hold on; title('Leaf On - Crown Polygons and associated points');
for j = 1:n_trees     
    % find raw las points within current polygon
    xpoly = laegernTreeTable_final.xPoly{j};
    ypoly = laegernTreeTable_final.yPoly{j};
    in = inpolygon(raw_on.x,raw_on.y,xpoly,ypoly);
    zpoly = raw_on.z(in);
    
    % calculate point cloud statistics 
    stats_on.zMax(j,1) = max(zpoly);         % max height
    stats_on.zMedian(j,1) = median(zpoly);   % median height
    stats_on.zMean(j,1) = mean(zpoly);       % mean height
    stats_on.zStd(j,1) = std(zpoly);         % standard deviation
    
end


%% save stats and pc structs to file

save(['/Users/scholl/geo_uzh/output/laegern/raw_' year '_leafoff_filtered.mat'], 'raw_off');
save(['/Users/scholl/geo_uzh/output/laegern/raw_' year '_leafon_filtered.mat'], 'raw_on');

save(['/Users/scholl/geo_uzh/output/laegern/stats_' year '_leafoff.mat'], 'stats_off');
save(['/Users/scholl/geo_uzh/output/laegern/stats_' year '_leafon.mat'], 'stats_on');


%% box plots of differences  leaf on - leaf off

stats_off_2010 = ;
stats_on_2010 = ; 
stats_off_2014 = ; 
stats_on_2014 = ; 

% restrict to the 8 species in Felix's paper



figure('Name',['Laegeren ' year ' Leaf Off - Leaf On Difference Boxplots']); 

subplot(4,1,1);
dif_zMedian =  stats_off.zMedian - stats_on.zMedian; 
boxplot(dif_zMedian,stats_off.species);title('median height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(4,1,2);
dif_zMean = stats_off.zMean - stats_on.zMean; 
boxplot(dif_zMean,stats_off.species);title('mean height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(4,1,3);
dif_zMax = stats_off.zMax - stats_on.zMax; 
boxplot(dif_zMax,stats_off.species);title('max height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(4,1,4);
dif_zStd = stats_off.zStd - stats_on.zStd; 
boxplot(dif_zStd,stats_off.species);title('std height'); 
line([0 20],[0 0],'color','k','linewidth',1)


