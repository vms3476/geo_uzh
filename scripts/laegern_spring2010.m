
% Spring 2010 Laegeren data

base = '/Users/scholl/geo_uzh/';
%dataPath = 'data/Laegeren/Spring_2010/PointCloud/';
%file = 'dat10041001009';
dataPath = 'data/Laegeren/Spring_2010/PointCloud/merged/';
file = 'merged_7_8_9_11';

cmap = parula;  % plot colormap
res = 1;      % raster resolution

% % read files using UZH methods
fileLAS = [base dataPath file '.las'];
% fileMAT = [base dataPath file '.mat'];
% pc_uzh = readlas(fileLAS); 
% mat = load(fileMAT); 
% data = mat.raw;

%% read using mparkan method
pc_mparkan = LASread(fileLAS);

% plot the pc
figure; 
myscatter3(pc_mparkan.record.x,pc_mparkan.record.y,pc_mparkan.record.z,pc_mparkan.record.z,cmap)

% keep only first returns 
i = ismember(pc_mparkan.record.return_number,1);
pc1.record.x = pc_mparkan.record.x(i);
pc1.record.y = pc_mparkan.record.y(i);
pc1.record.z = pc_mparkan.record.z(i);
pc1.record.intensity = pc_mparkan.record.intensity(i);
pc1.record.return_number = pc_mparkan.record.return_number(i);


% %% spatial subset based on Fabian's morphological traits map 
% outputFilepath = [base dataPath file '_subset.las'];
% xv = [669625; 669625; 670000; 670000]; 
% yv = [258900; 259270; 259270; 258900];
% s = LASclip(pc1, [xv, yv], outputFilepath, 'verbose', true);
% 
% % read the subset 
% pc_subset = LASread(outputFilepath); 

%% spatial subset based on Fabian's morphological traits map
x1 = 669625; 
x2 = 670000; 
y1 = 258900; 
y2 = 259270;
ii = pc1.record.x > (x1) & pc1.record.x < (x2) & pc1.record.y > (y1) & pc1.record.y < (y2);
pc1.record.xs = pc1.record.x(ii);
pc1.record.ys = pc1.record.y(ii);
pc1.record.zs = pc1.record.z(ii);

% plot 
figure; 
myscatter3(pc1.record.xs,pc1.record.ys,pc1.record.zs,pc1.record.zs,cmap)
title('Laegeren subset PC'); h = colorbar; xlabel(h,'Height [meters]','FontSize',11);


%% create raster of first returns 
ras1 = raw2ras([pc1.record.xs,pc1.record.ys,pc1.record.zs],0.5,0.5,'dsm'); 
[ras1.X, ras1.Y] = meshgrid(ras1.x,ras1.y);
figure; imagesc(ras1.x,ras1.y,ras1.z);

% read dtm 
load([base 'data/Laegeren/Spring_2010/DTM/dtm.mat']);

% % plot dtm
% figure; 
% myscatter3(dtm.X(:),dtm.Y(:),dtm.z(:),dtm.z(:), cmap);
% title('Laegeren dtm'); h = colorbar; xlabel(h,'meters','FontSize',11);

[dtm.X, dtm.Y] = meshgrid(dtm.x,dtm.y);
dtm.Z = interp2(dtm.X,dtm.Y,dtm.z,ras1.X,ras1.Y);
rasNorm = ras1.z - dtm.Z;

% plot normalized raster values
figure; 
myscatter3(ras1.X(:),ras1.Y(:),rasNorm(:),rasNorm(:),cmap)
title('Normalized raster - Laegern'); h = colorbar; xlabel(h,'Height [meters]','FontSize',11);
axis equal;axis tight;axis xy;
swisstick();
caxis([0,50])

figure; 
myscatter3(ras1.X(:),ras1.Y(:),rasNorm(:),rasNorm(:),cmap)
title('Normalized raster - Laegern'); h = colorbar; xlabel(h,'Height [meters]','FontSize',11);

hold on;

% plot polygons based on Fabian's tree morphology map.
% broad_vs_needle_laegern.m
for i = 1:numel(data_idx_deciduous)
    x = laegernTreeTable_final.xPoly{data_idx_deciduous(i)};
    y = laegernTreeTable_final.yPoly{data_idx_deciduous(i)};
    z = repmat(45,size(x));
    d = patch(x,y,z,[0 0.4 0.8]);
    d.FaceAlpha = 0.4;
end


for i = 1:numel(data_idx_coniferous)
    x = laegernTreeTable_final.xPoly{data_idx_coniferous(i)};
    y = laegernTreeTable_final.yPoly{data_idx_coniferous(i)};
        z = repmat(45,size(x));
    c = patch(x,y,z,[0, 0.8, 0.4]);
    c.FaceAlpha = 0.4;
end



%% write all las file extents in directory to shp file
outputFilepath = [base dataPath 'las_extents.shp'];
extent = LASextent([base dataPath],outputFilepath,'fig', true, 'verbose', true);

%% merge multiple las files

% it appears that tiles 7,8,9, and 10 will be useful for studying the
% region of Fabian's map

dataPath = 'data/Laegeren/Spring_2010/PointCloud/';

inputPath = {[base dataPath 'dat10041001007.las'], ...
             [base dataPath 'dat10041001008.las'], ...
             [base dataPath 'dat10041001009.las'], ...
             [base dataPath 'dat10041001011.las']};
%              [base dataPath 'dat10041001010.las'], ...


outputFilepath = [base dataPath 'merged/merged_7_8_9_11.las'];
s = LASmerge(inputPath, outputFilepath);

%% clip the raster based on the polygons. assess height difference between ? 



%% 