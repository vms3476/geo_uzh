
% specify main path to data
%base = '/media/Data/victoria/geo_uzh/';
base = '/Users/scholl/geo_uzh/';

% plot parameters
cmap = parula;% plot parameters

%% reading files 

path_data = 'data/KantonAargau/LeafOff/';
file = '624000_264000';

fname_las = [base path_data file '.laz'];
fname_mat = [base path_data file '_leaf_off.mat'];
fname_dtm = [base path_data 'DTM_' file '.tif']; 

% read mat file 
mat = load(fname_mat);
data = mat.data;

% % read las file
%raw = readlas(fname_las);
% % read DTM geotiff
% dtm = geotiffread(fname_dtm);

% to measure time elapsed
%tic
%elapsed = toc

%% subset tile

tic 

% find boundaries
x_min = min(data.x);
x_max = max(data.x);
y_min = min(data.y);
y_max = max(data.y);

% find indices based on specified boundaries
i = data.x > (x_min + 600) & data.x < (x_min + 700) & data.y > (y_min + 400) & data.y < (y_min + 500);
subset.x = data.x(i);
subset.y = data.y(i);
subset.z = data.z(i);
subset.z_AG = data.z_AG(i);
subset.rnnr = data.rnnr(i);
subset.z_DTM = data.z_DTM(i); 

% DSM using all points
subset_ras = raw2ras([subset.x,subset.y,subset.z],0.5,0.5,'dsm');
%save([root 'output/' 'dsm_subset_' file],'subset_ras')
[subset_ras.X, subset_ras.Y] = meshgrid(subset_ras.x,subset_ras.y);

% first return indices
r = subset.rnnr;
ii = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;

% DSM using only first returns
subset.x_r1 = subset.x(ii);
subset.y_r1 = subset.y(ii);
subset.z_r1 = subset.z(ii);
subset_r1_ras = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],0.5,0.5,'dsm');
%save([root path_data 'dsm_subset_r1_' file],'subset_r1_ras')
[subset_r1_ras.X, subset_r1_ras.Y] = meshgrid(subset_r1_ras.x,subset_r1_ras.y);

% subtract DTM from DSM of first returns 
subset_DTM_ras = raw2ras([subset.x,subset.y,subset.z_DTM],0.5,0.5,'dsm');
subset_dif_r1_dtm = subset_r1_ras.z - subset_DTM_ras.z;

% subset first returns from z_AG
subset.z_AG_r1 = subset.z_AG(ii);

elapsed_time = toc

%% plot 

% all points within subset
figure;  
myscatter3(subset.x, subset.y, subset.z, subset.z, cmap);
title('z values'); swisstick

%  z_AG
figure; 
myscatter3(subset.x, subset.y, subset.z_AG, subset.z_AG, cmap);
title('z AG values');swisstick 

% DTM 
figure; 
myscatter3(subset.x, subset.y, subset.z_DTM, subset.z_DTM, cmap);
title('z DTM values'); swisstick

% DSM (all returns)
figure; 
myscatter3(subset_ras.X(:), subset_ras.Y(:), subset_ras.z(:), subset_ras.z(:), cmap);
title('DSM, all points'); swisstick

% DSM (first returns)
figure;
myscatter3(subset_r1_ras.X(:), subset_r1_ras.Y(:), subset_r1_ras.z(:), subset_ras.z(:), cmap);
title('DSM, first returns'); swisstick

% DSM (first returns) - DTM
figure;
myscatter3(subset_r1_ras.X(:), subset_r1_ras.Y(:),subset_dif_r1_dtm(:),subset_dif_r1_dtm(:),cmap);
title('DSM - DTM, first returns only'); swisstick

% z AG (first returns)
figure;
myscatter3(subset.x_r1(:), subset.y_r1(:),subset.z_AG_r1(:),subset.z_AG_r1(:),cmap);
title('z AG, first returns only'); swisstick

%% plot everything in a single plot

figure; 

% all points within subset
subplottight(2,3,1);           
myscatter3(subset.x, subset.y, subset.z, subset.z, cmap);
title('z values'); swisstick

% z_AG
% subplottight(2,3,2);            
% myscatter3(subset.x, subset.y, subset.z_AG, subset.z_AG, cmap);
% title('z AG values');swisstick 

% DTM
subplottight(2,3,3);
myscatter3(subset.x, subset.y, subset.z_DTM, subset.z_DTM, cmap);
title('z DTM values'); swisstick

% DSM using all points
subplottight(2,3,4);
myscatter3(subset_ras.X(:), subset_ras.Y(:), subset_ras.z(:), subset_ras.z(:), cmap);
title('DSM, all points'); swisstick

% DSM using first returns
subplottight(2,3,5);
myscatter3(subset_r1_ras.X(:), subset_r1_ras.Y(:), subset_r1_ras.z(:), subset_ras.z(:), cmap);
title('DSM, first returns'); swisstick

% DSM (first returns) - DTM
subplottight(2,3,6);
myscatter3(subset_r1_ras.X(:), subset_r1_ras.Y(:),subset_dif_r1_dtm(:),subset_dif_r1_dtm(:),cmap);
title('DSM - DTM, first returns only'); swisstick

% z AG (first returns)
subplottight(2,3,2);
myscatter3(subset.x_r1(:), subset.y_r1(:),subset.z_AG_r1(:),subset.z_AG_r1(:),cmap);
title('z AG, first returns only'); swisstick

%% raster resolution

% 0.5
subset_r1_ras_05 = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],0.5,0.5,'dsm');
[subset_r1_ras_05.X, subset_r1_ras_05.Y] = meshgrid(subset_r1_ras_05.x,subset_r1_ras_05.y);

% 1.0v
subset_r1_ras_1 = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],1.0,1.0,'dsm');
[subset_r1_ras_1.X, subset_r1_ras_1.Y] = meshgrid(subset_r1_ras_1.x,subset_r1_ras_1.y);

% 2.0
subset_r1_ras_2 = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],2.0,2.0,'dsm');
[subset_r1_ras_2.X, subset_r1_ras_2.Y] = meshgrid(subset_r1_ras_2.x,subset_r1_ras_2.y);

figure; 
subplottight(1,3,1);
myscatter3(subset_r1_ras_05.X(:), subset_r1_ras_05.Y(:),subset_r1_ras_05.z(:),subset_r1_ras_05.z(:),cmap);
title('raster res: 0.5');

subplottight(1,3,2);
myscatter3(subset_r1_ras_1.X(:), subset_r1_ras_1.Y(:),subset_r1_ras_1.z(:),subset_r1_ras_1.z(:),cmap);
title('raster res: 01.0');

subplottight(1,3,3);
myscatter3(subset_r1_ras_2.X(:), subset_r1_ras_2.Y(:),subset_r1_ras_2.z(:),subset_r1_ras_2.z(:),cmap);
title('raster res: 2.0');

%% MATLAB coordinate test

file = '658000_245000';
leafoff = load([base 'data/KantonAargau/LeafOff/658000_245000_leaf_off.mat']);
leafon = load([base 'data/KantonAargau/LeafOn/658000_245000_leaf_on.mat']);

% find boundaries
x_min = min(leafoff.data.x);
x_max = max(leafoff.data.x);
y_min = min(leafoff.data.y);
y_max = max(leafoff.data.y);

