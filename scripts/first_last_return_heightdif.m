
%% load laegern 2010 leaf off data and dtm over aoi
raw_off_all = readlas('/Users/scholl/geo_uzh/data/Laegeren/2014_off/raw_2014_off.las');
dtmFilepath = '/Users/scholl/geo_uzh/data/Laegeren/2010_off/DTM/dtm.mat';
load(dtmFilepath); 
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y); 
raw_off_all.tz = raw_off_all.z;
% normalize PC to DTM
raw_off_all.z = raw_off_all.tz - interp2(dtm.X,dtm.Y,dtm.z,raw_off_all.x,raw_off_all.y);
% plot 
figure; myscatter3(raw_off_all.x,raw_off_all.y,raw_off_all.z,raw_off_all.z,parula); hold on;
title(['2010 Leaf Off PC Normalized']); colorbar; %swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on
load('/Users/scholl/geo_uzh/data/Fabian/laegernTreeTable_final20160629.mat')
overlay_polygons(laegernTreeTable_final)



%%

% subset that contains deciduous and coniferous trees
xMin = 6.698246099000014e+05; xMax = 6.698622806000002e+05;
yMin = 2.590395454000011e+05; yMax = 2.590693898999989e+05;
j = raw_off_all.x > (xMin) & raw_off_all.x < (xMax) & raw_off_all.y > (yMin) & raw_off_all.y < (yMax);
x = raw_off_all.x(j);
y = raw_off_all.y(j);
z = raw_off_all.z(j);
rnnr = raw_off_all.rnnr(j)';
figure; myscatter3(x,y,z,z,parula);

% entire tile 
x = raw_off_all.x';
y = raw_off_all.y';
z = raw_off_all.z';
rnnr = raw_off_all.rnnr';
figure; myscatter3(x,y,z,z,parula);

% remove middle returns
i = rnnr == 32 | rnnr == 42 | rnnr == 43 | rnnr == 52 | rnnr == 53 | rnnr == 54 | ...
    rnnr == 62 | rnnr == 63 | rnnr == 64 | rnnr == 65 | rnnr == 72 | rnnr == 73 | ...
    rnnr == 74 | rnnr == 75 | rnnr == 76 | rnnr == 11; 

x(i) = [];
y(i) = [];
z(i) = [];
rnnr(i) = [];

drnnr = diff(rnnr);
dz = diff(z);

% find indices of first returns for first-last echo pairs
% only keep negative differences
ii = (drnnr == 1 | drnnr == 2 | drnnr == 3 | drnnr == 4 | drnnr == 5 | drnnr == 6) & dz<0; 
xDif = x(ii);
yDif = y(ii); 
zDif = -dz(ii); 

% plot first-last return data with crown polygons overlaid
figure; myscatter3(xDif,yDif,zDif,zDif,parula); swisstick;
overlay_polygons(laegernTreeTable_final)
title(['2010 Leaf Off Laegeren'],'FontSize',14); c = colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 round(max(zDif))]); grid on;
ylabel(c,'First-Last Echo Heigh Difference [m]','FontSize',14)

% for testing, vertically stack the return number, diff(return number), 
% indices of first return pairs to keep for analysis, difference values
a = zeros(4,size(x,2),'double'); 
a(1,:) = rnnr;
a(2,:) = [drnnr 0.0];
a(3,:) = [ii 0.0];
a(4,:) = [dz 0.0];

% determine average first-last pulse distance difference per polygon
n_trees = size(laegernTreeTable_final,1);
stats.idField = laegernTreeTable_final.idField;
stats.species = laegernTreeTable_final.species;
for j = 1:n_trees     
    
    % find raw las points within current polygon
    xpoly = laegernTreeTable_final.xPoly{j};
    ypoly = laegernTreeTable_final.yPoly{j};
    in = inpolygon(xDif,yDif,xpoly,ypoly);
    zpoly = zDif(in);
    
    stats.zDif(j,1) = mean(zpoly); 
end

% keep only species 11 14 22 23 29 31 56 59, each has 40 or more polygons
k = ismember(stats.species,[11 14 22 23 29 31 56 59]); 
stats_plot = subsetraw(stats,k);
figure; boxplot(stats_plot.zDif,stats_plot.species);title('2010 Laegeren first-last echo height difference per pulse','FontSize',14); 
xlabel('Species','FontSize',14); ylabel('average height difference within crown polygon','FontSize',14);