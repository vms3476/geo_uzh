%% Testing leaf on/off, deciduous/coniferous

% specify main path to data
base = '/Users/scholl/geo_uzh/';

% plot parameters
cmap = parula;% plot parameters

% raster resolution
%res = 0.5;
%res = 1.0; 
res = 2.0; 

% read data

% leaf off
filename = [base 'data/KantonAargau/LeafOff/631000_264000_leaf_off.mat'];
out = filename(end-25:end-4);
tile = filename(end-25:end-13);
mat_leafoff = load(filename);
data_leafoff = mat_leafoff.data;
pc_leafoff = readlas([base 'data/KantonAargau/LeafOff/631000_264000.laz']);

% define subset boundaries for deciduous and coniferous

% % default: minimum x,y coordinates (lower left corner)
% x_min = min(data.x);
% x_max = max(data.x);
% y_min = min(data.y);
% y_max = max(data.y);

% deciduous
tree_type = 'deciduous';
x_min = 631760;
x_max = 631860;
y_min = 264610;
y_max = 264710;


% create and save mat file of dsm-dtm 

% define subset by x,y min,max from previous section
i = data_leafoff.x > (x_min) & data_leafoff.x < (x_max) & data_leafoff.y > (y_min) & data_leafoff.y < (y_max);

subset.x = data_leafoff.x(i);
subset.y = data_leafoff.y(i);
subset.z = data_leafoff.z(i);
subset.z_AG = data_leafoff.z_AG(i);
subset.rnnr = data_leafoff.rnnr(i);
subset.z_DTM = data_leafoff.z_DTM(i); 

% first return indices
r = subset.rnnr; ii = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;

% DSM using only first returns
subset.x_r1 = subset.x(ii);
subset.y_r1 = subset.y(ii);
subset.z_r1 = subset.z(ii);
subset_r1_ras = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],res,res,'dsm');
[subset_r1_ras.X, subset_r1_ras.Y] = meshgrid(subset_r1_ras.x,subset_r1_ras.y);

% subtract DTM from DSM of first returns 
dtm = raw2ras([subset.x,subset.y,subset.z_DTM],res,res,'dsm');
subset_r1_ras.dif = subset_r1_ras.z - dtm.z;

save([base 'output/' out '_' tree_type '_subset_res_' strrep(num2str(res),'.','')],'subset_r1_ras')


% coniferous 
tree_type = 'coniferous';
x_min = 631840;
x_max = 631940;
y_min = 264300;
y_max = 264400;

% create and save mat file of dsm-dtm 

% define subset by x,y min,max from previous section
i = data_leafoff.x > (x_min) & data_leafoff.x < (x_max) & data_leafoff.y > (y_min) & data_leafoff.y < (y_max);

subset.x = data_leafoff.x(i);
subset.y = data_leafoff.y(i);
subset.z = data_leafoff.z(i);
subset.z_AG = data_leafoff.z_AG(i);
subset.rnnr = data_leafoff.rnnr(i);
subset.z_DTM = data_leafoff.z_DTM(i); 

% first return indices
r = subset.rnnr; ii = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;

% DSM using only first returns
subset.x_r1 = subset.x(ii);
subset.y_r1 = subset.y(ii);
subset.z_r1 = subset.z(ii);
subset_r1_ras = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],res,res,'dsm');
[subset_r1_ras.X, subset_r1_ras.Y] = meshgrid(subset_r1_ras.x,subset_r1_ras.y);

% subtract DTM from DSM of first returns 
dtm = raw2ras([subset.x,subset.y,subset.z_DTM],res,res,'dsm');
subset_r1_ras.dif = subset_r1_ras.z - dtm.z;

save([base 'output/' out '_' tree_type '_subset_res_' strrep(num2str(res),'.','')],'subset_r1_ras')


%% leaf on 
filename = [base 'data/KantonAargau/LeafOn/631000_264000_leaf_on.mat'];
out = filename(end-24:end-4);
tile = filename(end-24:end-12);
mat_leafon = load(filename);
data_leafon = mat_leafon.data;
pc_leafon = readlas([base 'data/KantonAargau/LeafOn/631000_264000.laz']);

% deciduous
tree_type = 'deciduous';
x_min = 631760;
x_max = 631860;
y_min = 264610;
y_max = 264710;

% create and save mat file of dsm-dtm 

% define subset by x,y min,max from previous section
i = data_leafon.x > (x_min) & data_leafon.x < (x_max) & data_leafon.y > (y_min) & data_leafon.y < (y_max);

subset.x = data_leafon.x(i);
subset.y = data_leafon.y(i);
subset.z = data_leafon.z(i);
subset.z_AG = data_leafon.z_AG(i);
subset.rnnr = data_leafon.rnnr(i);
subset.z_DTM = data_leafon.z_DTM(i); 

% first return indices
r = subset.rnnr; ii = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;

% DSM using only first returns
subset.x_r1 = subset.x(ii);
subset.y_r1 = subset.y(ii);
subset.z_r1 = subset.z(ii);
subset_r1_ras = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],res,res,'dsm');
[subset_r1_ras.X, subset_r1_ras.Y] = meshgrid(subset_r1_ras.x,subset_r1_ras.y);

% subtract DTM from DSM of first returns 
dtm = raw2ras([subset.x,subset.y,subset.z_DTM],res,res,'dsm');
subset_r1_ras.dif = subset_r1_ras.z - dtm.z;

save([base 'output/' out '_' tree_type '_subset_res_' strrep(num2str(res),'.','')],'subset_r1_ras')


% coniferous 
tree_type = 'coniferous';
x_min = 631840;
x_max = 631940;
y_min = 264300;
y_max = 264400;

% create and save mat file of dsm-dtm 

% define subset by x,y min,max from previous section
i = data_leafon.x > (x_min) & data_leafon.x < (x_max) & data_leafon.y > (y_min) & data_leafon.y < (y_max);

subset.x = data_leafon.x(i);
subset.y = data_leafon.y(i);
subset.z = data_leafon.z(i);
subset.z_AG = data_leafon.z_AG(i);
subset.rnnr = data_leafon.rnnr(i);
subset.z_DTM = data_leafon.z_DTM(i); 

% first return indices
r = subset.rnnr; ii = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;

% DSM using only first returns
subset.x_r1 = subset.x(ii);
subset.y_r1 = subset.y(ii);
subset.z_r1 = subset.z(ii);
subset_r1_ras = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],res,res,'dsm');
[subset_r1_ras.X, subset_r1_ras.Y] = meshgrid(subset_r1_ras.x,subset_r1_ras.y);

% subtract DTM from DSM of first returns 
dtm = raw2ras([subset.x,subset.y,subset.z_DTM],res,res,'dsm');
subset_r1_ras.dif = subset_r1_ras.z - dtm.z;

save([base 'output/' out '_' tree_type '_subset_res_' strrep(num2str(res),'.','')],'subset_r1_ras')




%% plot both leaf off / on 

% load leaf on / off, deciduous / coniferous for the current tile

load([base 'output/' tile '_leaf_on_deciduous_subset_res_' strrep(num2str(res),'.','') '.mat']); leafon_deciduous = subset_r1_ras;
load([base 'output/' tile '_leaf_on_coniferous_subset_res_' strrep(num2str(res),'.','') '.mat']); leafon_coniferous = subset_r1_ras;
load([base 'output/' tile '_leaf_off_deciduous_subset_res_' strrep(num2str(res),'.','') '.mat']); leafoff_deciduous = subset_r1_ras;
load([base 'output/' tile '_leaf_off_coniferous_subset_res_' strrep(num2str(res),'.','') '.mat']); leafoff_coniferous = subset_r1_ras;


%z_max = round(max(max(max(leafoff.dif)),max(max(leafon.dif))));

% deciduous leaf on vs. leaf off DSM 
figure('Name',['Deciduous ' tile ' DSM - DTM, first returns only']); 
subplot(1,2,1);
myscatter3(leafoff_deciduous.X(:), leafoff_deciduous.Y(:),leafoff_deciduous.dif(:),leafoff_deciduous.dif(:),cmap);
title('Leaf Off'); swisstick
% zlim([0,z_max]);
% colorbar

subplot(1,2,2);
myscatter3(leafon_deciduous.X(:), leafon_deciduous.Y(:),leafon_deciduous.dif(:),leafon_deciduous.dif(:),cmap);
title('Leaf On'); swisstick
% zlim([0,z_max]);
% colorbar


% coniferous leaf on vs. leaf off DSM 
figure('Name',['Coniferous ' tile ' DSM - DTM, first returns only']); 
subplot(1,2,1);
myscatter3(leafoff_coniferous.X(:), leafoff_coniferous.Y(:),leafoff_coniferous.dif(:),leafoff_coniferous.dif(:),cmap);
title('Leaf Off '); swisstick
% zlim([0,z_max]);
% colorbar

subplot(1,2,2);
myscatter3(leafon_coniferous.X(:), leafon_coniferous.Y(:),leafon_coniferous.dif(:),leafon_coniferous.dif(:),cmap);
title('Leaf On'); swisstick
% zlim([0,z_max]);
% colorbar

