% create raster using leaf-off first returns for broad vs. needle leaf 

base = '/Users/scholl/geo_uzh/';
dataPath = 'data/KantonAargau/LeafOff/';
file = '631000_264000';

cmap = parula;  % plot colormap
res = 1;      % raster resolution

% read files
fileLAS = [base dataPath file '.laz'];
fileMAT = [base dataPath file '_leaf_off.mat'];
fileDTM = [base dataPath 'DTM_' file '.tif']; 

pc = readlas(fileLAS); 
mat = load(fileMAT); data = mat.data;
dtm = geotiffread(fileDTM);



%% subset forest types

% broad leaf

x1 = 631760;
x2 = 631860;
y1 = 264610;
y2 = 264710;
ii = data.x > (x1) & data.x < (x2) & data.y > (y1) & data.y < (y2);
broadleaf.x =  data.x(ii);
broadleaf.y =  data.y(ii);
broadleaf.z =  data.z(ii);
broadleaf.rnnr =  data.rnnr(ii);
broadleaf.z_DTM = data.z_DTM(ii);

% coniferous %tree_type = 'coniferous';
x1 = 631840;
x2 = 631940;
y1 = 264300;
y2 = 264400;
ii = data.x > (x1) & data.x < (x2) & data.y > (y1) & data.y < (y2);
evergreen.x =  data.x(ii);
evergreen.y =  data.y(ii);
evergreen.z =  data.z(ii);
evergreen.rnnr =  data.rnnr(ii);
evergreen.z_DTM = data.z_DTM(ii);


% normalize point clouds
broadleaf.z_norm = broadleaf.z - broadleaf.z_DTM;
evergreen.z_norm = evergreen.z - evergreen.z_DTM;


%% keep only first echos

r = broadleaf.rnnr;
i = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;
broadleaf.x1 = broadleaf.x(i);
broadleaf.y1 = broadleaf.y(i);
broadleaf.z1 = broadleaf.z(i); 
broadleaf.rnnr1 = broadleaf.rnnr(i); 
broadleaf.z_DTM1 = broadleaf.z_DTM(i); 
broadleaf.z_norm1 = broadleaf.z_norm(i);


r = evergreen.rnnr;
i = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;
evergreen.x1 = evergreen.x(i);
evergreen.y1 = evergreen.y(i);
evergreen.z1 = evergreen.z(i); 
evergreen.rnnr1 = evergreen.rnnr(i); 
evergreen.z_DTM1 = evergreen.z_DTM(i); 
evergreen.z_norm1 = evergreen.z_norm(i);


%% rasterize first echos 
rasBroadleaf = raw2ras([broadleaf.x1,broadleaf.y1,broadleaf.z_norm1],res,res,'dsm');
[rasBroadleaf.X,rasBroadleaf.Y] = meshgrid(rasBroadleaf.x,rasBroadleaf.y);

rasEvergreen = raw2ras([evergreen.x1,evergreen.y1,evergreen.z_norm1],res,res,'dsm');
[rasEvergreen.X,rasEvergreen.Y] = meshgrid(rasEvergreen.x,rasEvergreen.y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% plot 


% tile subset rasters 
figure; 
myscatter3(rasBroadleaf.X(:),rasBroadleaf.Y(:),rasBroadleaf.z(:), rasBroadleaf.z(:), cmap);
title('Broadleaf Subset - First Echos - Raster'); swisstick; h = colorbar; xlabel(h,'meters','FontSize',11);


figure; 
myscatter3(rasEvergreen.X(:),rasEvergreen.Y(:),rasEvergreen.z(:), rasBroadleaf.z(:), cmap);
title('Evergreen Subset - First Echos - Raster'); swisstick; h = colorbar; xlabel(h,'meters','FontSize',11);


% tile subset point cloud first echos 
figure; 
myscatter3(broadleaf.x1(:), broadleaf.y1(:), broadleaf.z_norm1(:), broadleaf.z_norm1(:), cmap);
title('Broadleaf Subset - First Echos - PC'); swisstick; h = colorbar; xlabel(h,'meters','FontSize',11);

figure; 
myscatter3(evergreen.x1(:), evergreen.y1(:), evergreen.z_norm1(:), evergreen.z_norm1(:), cmap);
title('Evergreen Subset - First Echos - PC'); swisstick; h = colorbar; xlabel(h,'meters','FontSize',11);





% % entire tile, all echos
% figure; 
% scatter3(data.x(:), data.y(:), data.z(:), 6, data.z(:), 'Marker', '.');
% colorbar; axis equal tight vis3d
% title('Tile first echos'); xlabel('X'); ylabel('Y'); ylabel('Z')


%% tile subsets, all echos 

figure; 
myscatter3(broadleaf.x(:), broadleaf.y(:), broadleaf.z_norm(:), broadleaf.z_norm(:),cmap);
title('Broadleaf Subset - All Echos - PC'); swisstick; h = colorbar; xlabel(h,'meters','FontSize',11);


figure; 
myscatter3(evergreen.x(:), evergreen.y(:), evergreen.z_norm(:), evergreen.z_norm(:), cmap);
title('Evergreen Subset - All Echos - PC'); swisstick; h = colorbar; xlabel(h,'meters','FontSize',11);




% %% Leaf on for comparison 
% base = '/Users/scholl/geo_uzh/';
% dataPath = 'data/KantonAargau/LeafOn/';
% fileOn = '631000_264000';
% 
% cmap = parula;  % plot colormap
% res = 1;      % raster resolution
% 
% % read files
% fileLASOn = [base dataPath file '.laz'];
% fileMATOn = [base dataPath file '_leaf_on.mat'];
% fileDTMOn = [base dataPath 'DTM_' file '.tif']; 
% 
% pcOn = readlas(fileLASOn); 
% matOn = load(fileMATOn); dataOn = matOn.data;
% dtmOn = geotiffread(fileDTMOn);

%% Laegeren 

% DSM first returns 
load([base 'data/Laegeren/Spring_2010/DSM_First_Echo/DSM_2010_Leafoff.mat']);
figure; 
[dsm10.X,dsm10.Y] = meshgrid(dsm10.x,dsm10.y);
myscatter3(dsm10.X(:),dsm10.Y(:),dsm10.z(:),dsm10.z(:),cmap);


