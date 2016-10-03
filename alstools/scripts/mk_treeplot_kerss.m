cd /Users/morsdorf/Projects/Active/Hyperforest/Cornelia/Singletreeextraction 
close all
load /Users/morsdorf/Projects/Active/Hyperforest/Cornelia/Results_Fco_Lai/VoxelGrids/LAIvoxelGrid10m.mat
ii = isnan(paiMatrix3d);
pai = paiMatrix3d;
pai(ii) = 0;
ii = isinf(pai);
pai(ii) = 0;
x = xVector(1:148);
y = yVector(1:211);
z = zVector(1:14);
[X,Y,Z] = meshgrid(x,y,z);
ii = pai > 0;
myscatter3(X(ii),Y(ii),Z(ii),pai(ii),flipud(summer(256)));

return
ltree.x = [];
ltree.y = [];
ltree.h = [];
ltree.dia = [];
close all
for i = 1:4
load(['Quartal',num2str(i),'_tree.mat']);
trees.dia(trees.dia >10) = 10; 
trees.h(trees.h >45) = 45; 
scatter(trees.x,trees.y,trees.dia*4,trees.h);
ltree.x = [ltree.x,trees.x];
ltree.y = [ltree.y,trees.y];
ltree.h = [ltree.h,trees.h];
ltree.dia = [ltree.dia,trees.dia];
hold on
axis equal;axis tight;
swisstick
colormap(ocean2);
colorbar
end

length(ltree.x)
S = shaperead('/Users/morsdorf/Projects/Active/Hyperforest/Data/kers/Hyperforest_Kersselaerspleyn_KV+16Cirkels.shp');
hold on
ii = [S.DBH_mm] > 120;
ftree.x = [S(ii).X];
ftree.y = [S(ii).Y];
ftree.h = [S(ii).Height_m];
ftree.dbh = [S(ii).DBH_mm];
figure;
are = [min(ftree.x) max(ftree.x) min(ftree.y) max(ftree.y)];
ii = ltree.x >= are(1) & ltree.x <= are(2) & ltree.y >= are(3) & ltree.y <= are(4);
ltree = subsetraw(ltree,ii);
scatter(ftree.x,ftree.y,ftree.dbh/100,ftree.h);
hold on
axis equal;axis tight;
swisstick
colormap(ocean2);
colorbar

trees = matchtrees(ftree,ltree);
load ../StudyArea.mat
load /Users/morsdorf/Projects/Active/Hyperforest/kersselaerspleyn.mat
ii = raw.x >= are(1) & raw.x <= are(2) & raw.y >= are(3) & raw.y <= are(4);
figure;
raw = subsetraw(raw,ii);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
raw.tz = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
clf; myscatter3(raw.x,raw.y,raw.tz,raw.tz,flipud(summer(256)));
hold on;maxh = nanmax(raw.tz);
plot3([ftree.x;ftree.x],[ftree.y;ftree.y],[ones(size([ftree.x]))*0;ones(size([ftree.x]))*maxh],'-k','linewidth',4);
plot3([ftree.x],[ftree.y],[ones(size([ftree.x]))*0],'*r','markersize',10);
plot3([ftree.x],[ftree.y],[ones(size([ftree.x]))*maxh],'*r','markersize',10);
%mat3d2osg('pointcloud_kers_stems');


chm = raw2ras([raw.x,raw.y,raw.tz],0.5,0.5,'dsm');
chm.z = inpaint_nans(chm.z,4);
[ii,chms] = locmax(chm.z,6,1.5,3);
figure;
myimage(chm.x,chm.y,chms);
hold on
[chm.X,chm.Y] = meshgrid(chm.x,chm.y);
plot(chm.X(ii),chm.Y(ii),'.w');
plot([S.X],[S.Y],'k.');
loc = [chm.X(ii);chm.Y(ii);chm.z(ii)];
%SZ = interp2(chm.X,chm.Y,chm.z,[S.X_crown],[S.Y_crown]);
RAW.x = raw.x;RAW.y = raw.y;RAW.z = raw.tz;
%loc = [[S.X_crown];[S.Y_crown];SZ];
%ii = [S.SOZ] == 1;
%loc = loc(:,ii);
%ii = ~isnan(loc(3,:));
%loc = loc(:,ii);
minheight = 3;
asp = 5;
method = 'gmmclust';

ntrees = segtree(RAW,loc',minheight,asp,method);