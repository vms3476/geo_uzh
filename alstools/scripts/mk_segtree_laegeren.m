% script for segmenting the forest at Laegeren into turbid media elements
cc
load /Volumes/Data1/Laegeren/aoi_spring.mat
load /Volumes/Data1/Laegeren/Spring_2010/DTM/dtm.mat
load /Volumes/Data1/Laegeren/Summer_2010/dsm_fe/dsm.mat
%load /Volumes/Data1/Laegeren/Summer_2010/DSMDTM.mat

raws = aoi_spring;
%are = [min(raws(:,1)) max(raws(:,1)) min(raws(:,2)) max(raws(:,2))];
% area for small TLS plot
are = [669660 669710 259040 259090];
ii = find(dsm.x >= are(1) & dsm.x <= are(2));
jj = find(dsm.y >= are(3) & dsm.y <= are(4));
dsm.x = double(dsm.x(ii));
dsm.y = double(dsm.y(jj));
dsm.z = double(dsm.z(jj,ii));
[dsm.X,dsm.Y] = meshgrid(dsm.x,dsm.y);
[dtm.X,dtm.Y] = meshgrid(dtm.x,dtm.y);
dtz = interp2(dtm.X,dtm.Y,dtm.z,dsm.X,dsm.Y);
dtm = dsm;dtm.z = dtz;clear dtz;
chm = dsm; chm.z = double(dsm.z)-double(dtm.z);
ii = raws(:,1) >= are(1) & raws(:,1) <= are(2) & raws(:,2) >= are(3) & raws(:,2) <= are(4);
raw.x = raws(ii,1);raw.y = raws(ii,2);raw.z = raws(ii,3) - interp2(dtm.X,dtm.Y,dtm.z,raw.x,raw.y);;
trees = segtree2(chm,raw,1,12,3);
save /Volumes/Data1/Laegeren/aoi_spring_segtree.mat