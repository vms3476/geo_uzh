% script for producing the best parameters for loctree

load /Volumes/Data2/ofenpass/gesamtgebiet/segtree2013/tiles/813_171_clust.mat
chm02 = chm;dsm02 = dsm;
load /Volumes/Data1/Ofenpass2010/segtree2013/dtmdsm2010.mat

dtm10 = subsetmod(dtm,chm02);
dsm10 = subsetmod(dsm,chm02);
chm10 = dsm10;
chm10.z = dsm10.z - dtm10.z
clear dtm dsm XI YI 
save /Volumes/Data1/Ofenpass2010/segtree2013/chm10 chm10 dsm10 dtm10
[chm10.X,chm10.Y] = meshgrid(chm10.x,chm10.y);
ii = locmax(chm.z,3,0.65,2);
x02 = chm.X(ii);
y02 = chm.Y(ii);

ii = locmax(chm10.z,5,0.65,2);
x10 = chm10.X(ii);
y10 = chm10.Y(ii);

myimage(chm10);
colormap(vegetation2)
hold on
plot(x02,y02,'.w','markersize',2);
plot(x10,y10,'.r','markersize',2);

