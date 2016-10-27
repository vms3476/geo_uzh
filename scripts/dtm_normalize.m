baseDir = ['/Users/scholl/geo_uzh/data/Laegeren/crown_map_subset/laegern_'];

%% 2010 data from flight lines
year = '2010';
load([baseDir year '_dtm.mat']);
[dtm.X,dtm.Y] = meshgrid(dtm.x, dtm.y); 

% normalize PC using dtm for leaf on and leaf off conditions
for l = {'off', 'on'}
    leaf = char(l)

    raw = readlas([baseDir year '_leaf' leaf '.las']);

    raw.ch = raw.z - interp2(dtm.X,dtm.Y,dtm.z, raw.x, raw.y); 

    raw = subsetraw(raw, raw.ch < 50);

    figure; myscatter3(raw.x,raw.y,raw.ch,raw.ch,parula);
    title([year ' Leaf ' leaf 'PC Normalized < 50m']); colorbar; swisstick
    axis equal;axis tight;axis xy; caxis([0 50]); grid on

    save([baseDir year '_leaf' leaf '_pc_norm.mat'], 'raw')

end


%% 2014 data is in tiles 
year = '2014'; 
dtm = load([baseDir year '_dtm.mat']);
dtm = dtm.data; 
[dtm.X,dtm.Y] = meshgrid(dtm.dtm_x, dtm.dtm_y); 

for l = {'off', 'on'}
    leaf = char(l)
    raw = readlas([baseDir year '_leaf' leaf '.las']);
    raw.ch = raw.z - interp2(dtm.X,dtm.Y,dtm.dtm, raw.x, raw.y); 
    raw = subsetraw(raw, raw.ch < 50);
    figure; myscatter3(raw.x,raw.y,raw.ch,raw.ch,parula);
    title([year ' Leaf ' leaf 'PC Normalized < 50m']); colorbar; swisstick
    axis equal;axis tight;axis xy; caxis([0 50]); grid on
    
end

%% 

% 2014 tile 1 258000
year = '2014';

leaf = 'on'

raw = readlas([baseDir year '_leaf' leaf '.las']);   
%rawFile = LASread([baseDir year '_leaf' leaf '.las']); %mparkan method
%raw = rawFile.record;

dtm1 = load([baseDir year '_669000_258000_dtm.mat']);
dtm1 = dtm1.data; 
[dtm1.X,dtm1.Y] = meshgrid(dtm1.dtm_x, dtm1.dtm_y); 
raw.ch = raw.z - interp2(dtm1.X,dtm1.Y,dtm1.dtm, raw.x, raw.y);
raw1 = subsetraw(raw, raw.ch < 50 & ~isnan(raw.ch));
figure; hold on; myscatter3(raw1.x,raw1.y,raw1.ch,raw1.ch,parula);
   

dtm2 = load([baseDir year '_669000_259000_dtm.mat']);
dtm2 = dtm2.data; 
[dtm2.X,dtm2.Y] = meshgrid(dtm2.dtm_x, dtm2.dtm_y); 
raw.ch = raw.z - interp2(dtm2.X,dtm2.Y,dtm2.dtm, raw.x, raw.y);
raw2 = subsetraw(raw, raw.ch < 50 & ~isnan(raw.ch));
myscatter3(raw2.x,raw2.y,raw2.ch,raw2.ch,parula);

title([year ' Leaf ' leaf 'PC Normalized < 50m']); colorbar; swisstick
axis equal;axis tight;axis xy; caxis([0 50]); grid on

% merge the normalized las data for both tiles
raw1.record = raw1;
raw2.record = raw2; 
inputLas = {raw2, raw1};
raw_2014_off = LASmerge(inputLas)

% merge 
fields = fieldnames(raw1);
for i = 1:length(fields)
merged.(fields{i}) = [raw1.(fields{i});raw2.(fields{i})];
end

figure; hold on; myscatter3(merged.x,merged.y,merged.ch,merged.ch,parula);

save([baseDir year '_leaf' leaf '_pc_norm.mat'],'merged')

%% Read in crown polygons

disp('Reading crown polygons...');

run crown_polygons_laegeren.m

%% statistics 
n_trees = size(laegernTreeTable_final,1);
stats_off.idField = laegernTreeTable_final.idField;
stats_off.species = laegernTreeTable_final.species;
stats_off.xPoly = laegernTreeTable_final.xPoly;
stats_off.yPoly = laegernTreeTable_final.yPoly;

% 2010

raw_off = load('/Users/scholl/geo_uzh/data/Laegeren/crown_map_subset/laegern_2010_leafoff_pc_norm.mat');
raw_off = raw_off.raw;


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
    stats_off.groundEchoFraction(j,1) = sum(g1) / numel(zpoly_above3m);
end

% leaf on
raw_on = load('/Users/scholl/geo_uzh/data/Laegeren/crown_map_subset/laegern_2010_leafon_pc_norm.mat');
raw_on = raw_on.raw;

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

%% difference

k = ismember(stats_off.species,[11 14 22 23 29 31 56 59]); 
stats.off = subsetraw(stats_off,k); 
stats.on = subsetraw(stats_on,k); 
stats.off.species = stats.species(k);

%% boxplot


% 2010 
figure('Name','Laegeren 2010 Leaf Off - Leaf On Difference Boxplots'); 

subplot(6,1,1);
dif_zMax = stats.off.zMax - stats.on.zMax; 
boxplot(dif_zMax,stats.off.species);title('max height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-3,3])

subplot(6,1,2);
dif_zMedian =  stats.off.zMedian - stats.on.zMedian; 
boxplot(dif_zMedian,stats.off.species);title('median height'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-4,3])

subplot(6,1,3);
dif_zStd = stats.off.zStd - stats.on.zStd; 
boxplot(dif_zStd,stats.off.species);title('std height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,4);
dif_zMean = stats.off.zMean - stats.on.zMean; 
boxplot(dif_zMean,stats.off.species);title('mean height'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,5);
dif_singleEchoFraction = stats.off.singleEchoFraction - stats.on.singleEchoFraction; 
boxplot(dif_singleEchoFraction,stats.off.species);title('fraction of single echos'); 
line([0 20],[0 0],'color','k','linewidth',1)

subplot(6,1,6);
dif_groundEchoFraction = stats.off.groundEchoFraction - stats.on.groundEchoFraction; 
boxplot(dif_groundEchoFraction,stats.off.species);title('fraction of ground echos'); 
line([0 20],[0 0],'color','k','linewidth',1); ylim([-0.2,0.4]);