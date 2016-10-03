%% MATLAB coordinate test

% specify main path to data
root = '/media/Data/victoria/geo_uzh/';

% plot parameters
cmap = parula;% plot parameters

% input mat filename

% leaf off
filename = [root 'data/KantonAargau/LeafOff/631000_264000_leaf_off.mat'];
out = filename(end-25:end-4);
tile = filename(end-25:end-13);

% leaf on
filename = [root 'data/KantonAargau/LeafOn/631000_264000_leaf_on.mat'];
out = filename(end-24:end-4);
tile = filename(end-24:end-13);

mat = load(filename);
data = mat.data;

%% define subset boundaries

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

% coniferous 
tree_type = 'coniferous';
x_min = 631840;
x_max = 631940;
y_min = 264300;
y_max = 264400;

% raster resolution
% res = 0.5;
res = 1.0; 

%% create and save mat file of dsm-dtm 

% find indices based on specified boundaries

% % define subset by number of pixels here
%i = data.x > (x_min + 200) & data.x < (x_min + 400) & data.y > (y_min + 200) & data.y < (y_min + 400);

% define subset by x,y min,max from previous section
i = data.x > (x_min) & data.x < (x_max) & data.y > (y_min) & data.y < (y_max);

subset.x = data.x(i);
subset.y = data.y(i);
subset.z = data.z(i);
subset.z_AG = data.z_AG(i);
subset.rnnr = data.rnnr(i);
subset.z_DTM = data.z_DTM(i); 

% first return indices
r = subset.rnnr;
ii = r==11 | r==21 | r==31 | r==41 | r==51 |r==61 | r==71;

% DSM using only first returns
subset.x_r1 = subset.x(ii);
subset.y_r1 = subset.y(ii);
subset.z_r1 = subset.z(ii);
subset_r1_ras = raw2ras([subset.x_r1,subset.y_r1,subset.z_r1],res,res,'dsm');
[subset_r1_ras.X, subset_r1_ras.Y] = meshgrid(subset_r1_ras.x,subset_r1_ras.y);

% subtract DTM from DSM of first returns 
dtm = raw2ras([subset.x,subset.y,subset.z_DTM],res,res,'dsm');
subset_r1_ras.dif = subset_r1_ras.z - dtm.z;

save([root 'output/' out '_' tree_type '_subset_res_' ],'subset_r1_ras')

%% plot both leaf off / on 

% load leaf on / off, deciduous / coniferous for the current tile
if findstr(out, 'on')
    tree_type = 'deciduous';
    load([root 'output/' out '_' tree_type '_subset']); leafon_deciduous = subset_r1_ras;
    load([root 'output/' out(1:end-1) 'ff' '_' tree_type '_subset']); leafoff_deciduous = subset_r1_ras;
    
    tree_type = 'coniferous';
    load([root 'output/' out '_' tree_type '_subset']); leafon_coniferous = subset_r1_ras;
    load([root 'output/' out(1:end-1) 'ff' '_' tree_type '_subset']); leafoff_coniferous = subset_r1_ras;
else
    tree_type = 'deciduous';
    load([root 'output/' out '_' tree_type '_subset']); leafoff_deciduous = subset_r1_ras;
    load([root 'output/' out(1:end-2) 'n' '_' tree_type '_subset']); leafon_deciduous = subset_r1_ras;
    
    tree_type = 'coniferous';
    load([root 'output/' out '_' tree_type '_subset']); leafoff_coniferous = subset_r1_ras;
    load([root 'output/' out(1:end-2) 'n' '_' tree_type '_subset']); leafon_coniferous = subset_r1_ras;
end


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

