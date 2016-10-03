close all
cd /Volumes/Data1/Belgium/kersselaerspleyn/
load /Volumes/Data1/Belgium/kersselaerspleyn/Kersselaerspleyn_raw.mat
load /Volumes/Data1/Belgium/kersselaerspleyn/Kersselaerspleyn.mat
S = shaperead('/Volumes/Data1/Belgium/FieldData/Hyperforest_KersselaerspleynKV+16Cirkels.shp');
ii = [S(:).DBH_mm]> 150;
x = [S(ii).X];y = [S(ii).Y];
ii = chm.x >= min(x) & chm.x <= max(x);
jj = chm.y >= min(y) & chm.y <= max(y);
chm.x = chm.x(ii);
chm.y = chm.y(jj);
chm.z = chm.z(jj,ii);
ii = raw.x >= min(x) & raw.x <= max(x) & raw.y >= min(y) & raw.y <= max(y);
raw = subsetraw(raw,ii);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
raw.z = raw.z - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);
chm.z(isnan(chm.z)) = 0;
minheight = 3;
filtsize = 9;
% test tree location algorithm
[tree] = locmax(chm.z,minheight,filtsize);
% get starting positions
[gx,gy] = meshgrid(chm.x,chm.y);
sx = gx(tree);
sy = gy(tree);
imagesc(chm.x,chm.y,chm.z);
hold on
swisstick;axis xy;axis equal;axis tight
plot(x,y,'.w','markersize',7);
plot(sx,sy,'.k','markersize',7);
print -depsc2 -r300 treelocs.eps
[x,y,z,idx,RAW] = segtree2(chm,raw,3,9,4);
trees = crdata(x,y,z);

