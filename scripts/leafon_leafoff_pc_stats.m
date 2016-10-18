% compare flight lines between 2010 and 2014

%% 2010 Spring from lidarlab/data/Laegeren/Spring_2010/Pointcloud/
dir2010_off = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/Pointcloud';

%% 2010 Summer from lidarlab/data/Laegeren/Summer_2010/Pointcloud/
dir2010_on = '/Users/scholl/geo_uzh/data/Laegeren/2010_on/Pointcloud/';

%% Merge and read las files within study area

% specify coordinates of study area
aoi = [669660 669960 258910 259210];

% 2010 Spring leaf off
raw2010_off = getrawlas(dir2010_off,aoi);

% 2010 Summer leaf on
raw2010_on = getrawlas(dir2010_on,aoi);

% Load DTM 
dtmFilepath = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/Pointcloud/DTM/dtm.mat';
load(dtmFilepath);

% interpolate dtm to the raw point coordinates
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
raw2010_off.tz = raw2010_off.z;
raw2010_off.z = raw2010_off.tz - interp2(dtm.X,dtm.Y,dtm.z,raw2010_off.x,raw2010_off.y);

raw2010_on.tz = raw2010_on.z;
raw2010_on.z = raw2010_on.tz - interp2(dtm.X,dtm.Y,dtm.z,raw2010_on.x,raw2010_on.y);

% plot, remove noisy values above 50m
i = raw2010_off.z>50;
figure; myscatter3(raw2010_off.x(i),raw2010_off.y(i),raw2010_off.z(i),raw2010_off.z(i),parula);
title('2010 Leaf Off PC Normalized'); colorbar

figure; myscatter3(raw2010_on.x,raw2010_on.y,raw2010_on.z,raw2010_on.z,parula);
title('2010 Leaf On PC Normalized'); colorbar
