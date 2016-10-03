% script for plotting tree location on klosters chm
%load /Volumes/Data1/KlostersDavosSNP/klosters_dsmdtm.mat
trees = importdata('/Volumes/Data1/KlostersDavos/Produkte/trees.csv',',',1); 
chm = dsm;
chm.z = dsm.z-dtm.z;
chm.z = chm.z(1:2:end,1:2:end);
chm.x = chm.x(1:2:end);
chm.y = chm.y(1:2:end);
myimage(chm);
caxis([0 50])
hold on
hp = plot(trees.data(:,1),trees.data(:,2),'.','color','k');
set(hp,'markersize',1)
colormap([[1,1,1];ocean2])
title('CHM Kloster/Davos and tree loactions');
colorbar
print -dtiff -r3600 chm_treeloc_davos.tif

